"""Parse and transform Qubit assay output files for LabKey ingestion."""

import re
# Pylint can't see ThreadPoolExecutor through concurrent.futures' re-export from
# concurrent.futures.thread — the import works at runtime.
from concurrent.futures import ThreadPoolExecutor  # pylint: disable=no-name-in-module
from datetime import datetime
from typing import Tuple, Dict, Any, List

import pandas as pd

from lk_sample_ID_assignment import (
    assign_sample_ids, query_labkey_rest_api, execute_query, PROD_HOST,
)

_BATCH_ID_RE = re.compile(r'^[A-Za-z0-9]+$')

ASSAY_MAP = {
    'HS': 'Total dsDNA High Sensitivity (HS)',
    'BR': 'Total dsDNA Broad Range (BR)',
}

DEVICE_MAP = {
    'L2': '00091',
    'L3': '00113',
}

SAMPLE_TYPE_CODES = {
    'Hu': 'Human DNA',
    'Bk': 'Bulk DNA',
    'Mb': 'Microbial DNA',
}

CONTROL_TERMS = ['POS', 'NEG', 'TEBuffer', 'NA12878']


def _fetch_operator(batch_id: str, host: str) -> str:
    """Fetch QubitOperator from LabKey job templates."""
    templates_df = execute_query(schema="sampleManagement.jobtypes", host=host)
    job_templates = templates_df['name'].tolist()

    def fetch_template(job_template):
        sql = f"""
        SELECT
        jt.Name,
        jt.BatchID,
        CAST(jt.QubitOperator.DisplayName AS VARCHAR) AS QubitOperator
        FROM {job_template} as jt
        WHERE jt.BatchID = '{batch_id}'
        """
        result = query_labkey_rest_api(
            "/query-executeSql.api", host=host,
            schemaName="sampleManagement.jobtypes", sql=sql,
        )
        return pd.DataFrame(result.get("rows", []))

    with ThreadPoolExecutor(max_workers=10) as executor:
        dfs = list(executor.map(fetch_template, job_templates))

    workflows_df = (
        pd.concat(dfs, axis=0, ignore_index=True)
        if dfs
        else pd.DataFrame(columns=["Name", "BatchID", "QubitOperator"])
    )
    workflows_df.fillna('-', inplace=True)
    return workflows_df['QubitOperator'].values[0] if len(workflows_df) > 0 else '-'


def _extract_batch_and_asset(qubit: pd.DataFrame) -> Tuple[pd.DataFrame, str, str]:
    """Parse Batch ID + asset code from `Sample ID` (preferred) or `Tags` (fallback).

    Picks the first row in which the source column contains a `|`, so standard /
    control rows at the top of the CSV don't hijack the extraction.
    """
    def _first_piped_row(col: str):
        series = qubit[col].astype(str)
        mask = series.str.contains(r'\|', regex=True, na=False)
        if not mask.any():
            return None
        return series.loc[mask].iloc[0]

    piped_value = None
    source_col = None
    for candidate in ('Sample ID', 'Tags'):
        if candidate in qubit.columns:
            piped_value = _first_piped_row(candidate)
            if piped_value is not None:
                source_col = candidate
                break

    if piped_value is None:
        raise ValueError(
            "CSV must contain either 'Sample ID' or 'Tags' with pipe-delimited Batch ID|AssetCode"
        )

    split = qubit[source_col].astype(str).str.split('|', n=1, expand=True)
    qubit['Batch ID'] = split[0].str.lstrip('#')
    qubit['Qubit Asset ID'] = split[1] if split.shape[1] > 1 else None

    batch_id, _, asset_code = piped_value.partition('|')
    batch_id = batch_id.lstrip('#')
    if not _BATCH_ID_RE.match(batch_id):
        raise ValueError(f"Invalid Batch ID '{batch_id}': must be alphanumeric (no spaces)")

    return qubit, batch_id, asset_code


def _extract_rundate(qubit: pd.DataFrame) -> str:
    """Extract a filesystem-safe run date stamp from the Test Date column."""
    series = qubit['Test Date'].dropna().astype(str)
    if series.empty:
        return 'unknown'
    raw = series.iloc[0]
    for fmt in ('%d/%m/%Y %I:%M:%S %p', '%Y-%m-%d %H:%M:%S', '%Y-%m-%d'):
        try:
            return datetime.strptime(raw, fmt).strftime('%Y%m%d_%H-%M-%S')
        except ValueError:
            continue
    return 'unknown'


def qubit_processing(
    qubit_data_file: str, host: str = PROD_HOST,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, List[str]], Dict[str, Any]]:
    """Process a Qubit CSV and return (results, run_info, pass_fail, metadata)."""
    qubit = pd.read_csv(qubit_data_file)
    qubit.fillna({'Reagent Lot#': "MISSING|MISSING", 'Tags': "MISSING|MISSING"}, inplace=True)

    qubit[['Qubit Standards Lot #', 'Qubit Dye Lot #']] = (
        qubit['Reagent Lot#'].str.split(r'\|', expand=True) + ' (Verified)'
    )

    qubit, batch_id, _ = _extract_batch_and_asset(qubit)

    assay_series = qubit['Assay Name'].dropna().astype(str)
    assay_code = assay_series.iloc[0][-2:] if not assay_series.empty else ''
    qubit['Assay Name'] = (
        qubit['Assay Name'].astype(str).str[-2:].map(ASSAY_MAP).fillna('INVALID ASSAY')
    )
    qubit['Qubit Asset ID'] = qubit['Qubit Asset ID'].map(DEVICE_MAP).fillna('INVALID ASSET ID')

    qubit['Sample Number'] = pd.to_numeric(
        qubit['Sample Name'].astype(str).str.extract(r'S(\d+)', expand=False),
        errors='coerce',
    )
    qubit = qubit.dropna(subset=['Sample Number']).copy()
    qubit['Sample Number'] = qubit['Sample Number'].astype(int)

    rundate_str = _extract_rundate(qubit)

    with ThreadPoolExecutor(max_workers=2) as executor:
        future_ids = executor.submit(assign_sample_ids, batch_id, host)
        future_operator = executor.submit(_fetch_operator, batch_id, host)

    id_map, sample_types_present = future_ids.result()
    operator = future_operator.result()

    qubit = id_map.merge(qubit, how='inner', left_on='SampleNumber', right_on='Sample Number')
    qubit['Operator'] = operator

    if qubit.empty:
        raise ValueError(
            "Manifest and DNA were found, but no rows matched the Sample Numbers "
            "(S1, S2, ...) parsed from 'Sample Name' in the Qubit CSV."
        )

    qubit = qubit[['Name', 'Original Sample Conc.', 'Std 1 RFU', 'Std 2 RFU', 'Sample RFU',
                   'Batch ID', 'Operator', 'Test Date', 'Qubit Standards Lot #',
                   'Qubit Dye Lot #', 'Assay Name', 'Qubit Asset ID']]
    qubit.rename(columns={
        'Name': 'Sample ID',
        'Original Sample Conc.': 'Sample Concentration',
        'Std 1 RFU': 'Standard1RFU',
        'Std 2 RFU': 'Standard2RFU',
    }, inplace=True)

    results_df = qubit[
        ['Sample ID', 'Sample Concentration', 'Standard1RFU', 'Standard2RFU', 'Sample RFU']
    ].reset_index(drop=True)

    run_info_df = qubit[
        ['Batch ID', 'Operator', 'Test Date', 'Qubit Standards Lot #',
         'Qubit Dye Lot #', 'Assay Name', 'Qubit Asset ID']
    ].head(1).reset_index(drop=True)

    control_mask = qubit['Sample ID'].str.contains(
        '|'.join(CONTROL_TERMS), na=False, regex=True
    )
    samp = qubit[~control_mask].copy()
    samp['Sample Concentration'] = pd.to_numeric(samp['Sample Concentration'], errors='coerce')
    bad = samp[samp['Sample Concentration'].isna()]
    if not bad.empty:
        offenders = ', '.join(bad['Sample ID'].astype(str).tolist())
        raise ValueError(
            f"Non-numeric 'Original Sample Conc.' value(s) found for sample(s): {offenders}. "
            f"Expected a numeric concentration. Fix the CSV."
        )
    pass_fail = {
        'Passed Samples (Conc >= 5)':
            samp.loc[samp['Sample Concentration'] >= 5, 'Sample ID'].tolist(),
        'Failed Samples (Conc < 5)':
            samp.loc[samp['Sample Concentration'] < 5, 'Sample ID'].tolist(),
    }

    sample_types_used = '_'.join(sorted(
        code for code, present in sample_types_present.items() if present
    ))

    metadata = {
        "batch_id": batch_id,
        "assay": assay_code,
        "rundate": rundate_str,
        "filename_stem": f"{rundate_str}_{assay_code}_{batch_id}_{sample_types_used}",
    }

    return results_df, run_info_df, pass_fail, metadata

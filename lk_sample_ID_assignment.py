"""Shared LabKey sample-ID assignment used by Origin Sciences converters.

Resolves a Batch ID to a `(SampleNumber, Sample ID)` mapping by querying the
SOP-98 and SOP-64 Manifest Creation assays plus the Human/Microbial/Bulk DNA
sample types, then merging on the OriCol aliquot. Lifted from the production
qubit.py converter so future converters can share one implementation.
"""

from concurrent.futures import ThreadPoolExecutor  # pylint: disable=no-name-in-module
from typing import Tuple, Dict, Any

import pandas as pd
import requests

PROD_HOST = "originsciences.app.labkey.host"
_FOLDER = "home"
_APIKEY = "5830a18d79e74c1589bda8df43762f1a348d0ad85d243c075465df6329fa4b29"
_REQUEST_TIMEOUT = 30

_session = requests.Session()
_session.headers["Authorization"] = f"LABKEY apikey={_APIKEY}"


def query_labkey_rest_api(
    path_controller_action: str, host: str = PROD_HOST, **params: Any,
) -> Dict[str, Any]:
    """HTTPS GET against a LabKey REST API endpoint, returns parsed JSON."""
    response = _session.get(
        f"https://{host}/{_FOLDER}/{path_controller_action}",
        params={"containerPath": _FOLDER, **params},
        timeout=_REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    return response.json()


def execute_query(schema: str = "core", host: str = PROD_HOST) -> pd.DataFrame:
    """List queries in a LabKey schema (wraps query-GetQueries.api)."""
    response = _session.get(
        f"https://{host}/{_FOLDER}/query-GetQueries.api",
        params={"schemaName": schema},
        timeout=_REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    return pd.DataFrame(response.json().get("queries", []))


def _addsuffix(df: pd.DataFrame, col: str, group: str) -> pd.DataFrame:
    """Add _N suffix to duplicate values in a column, leaving unique values unchanged."""
    df[col] += (
        df.groupby(group).cumcount().add(1).astype(str).radd('_')
        .mask(df.groupby(group)[col].transform('count') == 1, '')
    )
    return df


def _fetch_manifest_sql(manifest_name: str, batch_id: str, host: str) -> pd.DataFrame:
    # Use a schema-scoped query (schemaName=assay.General.<manifest_name>, FROM Data).
    # The three-part path `assay.General."<name>".Data` is rejected by LabKey on this
    # instance, and `BatchID` is not a direct column on Data — it's only on Run.
    sql = (
        'SELECT ProcessPosition, RowId, SampleID,'
        ' CAST(SampleID.Name AS VARCHAR) AS Sample_ID,'
        ' CAST(Run.BatchID AS VARCHAR) AS Batch_ID'
        ' FROM Data'
        f" WHERE Run.BatchID = '{batch_id}'"
        ' ORDER BY RowId'
    )
    result = query_labkey_rest_api(
        "query-executeSql.api", host=host,
        schemaName=f"assay.General.{manifest_name}", sql=sql,
        columns="", includeHidden=True,
    )
    return pd.DataFrame(result.get("rows", []))


def _fetch_manifest(batch_id: str, host: str) -> pd.DataFrame:
    """Fetch manifest rows from SOP-98 and SOP-64 Manifest Creation assays combined."""
    with ThreadPoolExecutor(max_workers=2) as executor:
        dfs = list(executor.map(
            lambda name: _fetch_manifest_sql(name, batch_id, host),
            ["SOP-98 Manifest Creation", "SOP-64 Manifest Creation"],
        ))

    combined = pd.concat(dfs, axis=0, ignore_index=True) if dfs else pd.DataFrame()
    if combined.empty or 'ProcessPosition' not in combined.columns:
        return pd.DataFrame(columns=['Sample_ID', 'SampleNumber'])

    combined['SampleNumber'] = combined['ProcessPosition'].str.extract(r'(\d+)').astype(int)
    return combined[['Sample_ID', 'SampleNumber']]


def _fetch_dna_table(assay_name: str, batch_column: str, batch_id: str, host: str) -> pd.DataFrame:
    """Fetch DNA samples from a LabKey sample type, with the batch-id column configurable."""
    sql = (
        f'SELECT Name, ExtractedFrom, {batch_column} AS BatchID,'
        ' CAST(ExtractedFrom.Name AS VARCHAR) AS ExtractedFrom_1'
        f' FROM samples."{assay_name}"'
        f" WHERE {batch_column} = '{batch_id}'"
        " ORDER BY Stored ASC"
    )
    result = query_labkey_rest_api(
        "query-executeSql.api", host=host,
        schemaName="core", sql=sql, columns="", includeHidden=True,
    )
    df = pd.DataFrame(result.get("rows", []))
    if df.empty:
        return pd.DataFrame(columns=['ExtractedFrom_1', 'Name'])
    return df[['ExtractedFrom_1', 'Name']]


def _fetch_bulk_dna(batch_id: str, host: str) -> pd.DataFrame:
    """Fetch Bulk DNA samples — uses SecondaryParent instead of ExtractedFrom."""
    sql = (
        'SELECT Name, SecondaryParent, ExperimentID AS BatchID,'
        ' CAST(SecondaryParent.Name AS VARCHAR) AS ExtractedFrom_1'
        f' FROM samples."DNA"'
        f" WHERE ExperimentID = '{batch_id}'"
        ' ORDER BY Stored ASC'
    )
    result = query_labkey_rest_api(
        "query-executeSql.api", host=host,
        schemaName="core", sql=sql, columns="", includeHidden=True,
    )
    df = pd.DataFrame(result.get("rows", []))
    if df.empty:
        return pd.DataFrame(columns=['ExtractedFrom_1', 'Name'])
    return df[['ExtractedFrom_1', 'Name']]


def assign_sample_ids(
    batch_id: str, host: str = PROD_HOST,
) -> Tuple[pd.DataFrame, Dict[str, bool]]:
    """Resolve LabKey Sample IDs for a given Batch ID.

    Fetches the SOP-98 + SOP-64 manifest plus Human/Microbial/Bulk DNA tables
    concurrently, dedupes duplicate `ExtractedFrom_1` / `Sample_ID` values with
    `_N` suffixes, then merges manifest ↔ DNA on the OriCol aliquot.

    Returns:
        (id_map, sample_types_present)
          id_map: DataFrame with columns ['SampleNumber', 'Name']
          sample_types_present: {'Hu': bool, 'Mb': bool, 'Bk': bool}
    """
    with ThreadPoolExecutor(max_workers=8) as executor:
        future_manifest = executor.submit(_fetch_manifest, batch_id, host)
        future_human = executor.submit(_fetch_dna_table, "Human DNA", "BatchID", batch_id, host)
        future_microbial = executor.submit(
            _fetch_dna_table, "Microbial DNA", "DEXID", batch_id, host,
        )
        future_bulk = executor.submit(_fetch_bulk_dna, batch_id, host)

    manifest = future_manifest.result()
    human_df = future_human.result()
    microbial_df = future_microbial.result()
    bulk_df = future_bulk.result()
    dna = pd.concat([human_df, microbial_df, bulk_df], ignore_index=True)

    if manifest.empty:
        raise ValueError(
            f"No manifest rows found for Batch ID '{batch_id}' in SOP-98 or SOP-64 "
            f"Manifest Creation. Check that the manifest exists in LabKey for this batch."
        )
    if dna.empty:
        raise ValueError(
            f"No DNA samples found for Batch ID '{batch_id}' in Human DNA (BatchID), "
            f"Microbial DNA (DEXID), or Bulk DNA (ExperimentID) in LabKey."
        )

    dna = _addsuffix(dna, 'ExtractedFrom_1', 'ExtractedFrom_1')
    manifest = _addsuffix(manifest, 'Sample_ID', 'Sample_ID')

    merged = manifest.merge(dna, how='inner', left_on='Sample_ID', right_on='ExtractedFrom_1')
    id_map = merged[['SampleNumber', 'Name']]

    sample_types_present = {
        'Hu': not human_df.empty,
        'Mb': not microbial_df.empty,
        'Bk': not bulk_df.empty,
    }

    return id_map, sample_types_present

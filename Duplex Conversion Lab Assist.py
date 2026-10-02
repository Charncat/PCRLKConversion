# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 2026

@author: DanielMoutia

Combined Duplex Conversion for QuantStudio 1 (.xls) and LightCycler 480 (.txt) exports.
Merges Duplex Conversion Excel.py and Duplex Conversion Notepad.py.
"""

# Import all packages (remove any that are not needed)

import re
import math
# Pylint can't see ThreadPoolExecutor through concurrent.futures' re-export from
# concurrent.futures.thread — the import works at runtime.
from concurrent.futures import ThreadPoolExecutor  # pylint: disable=no-name-in-module

import pandas as pd

from backend.labkey_interface import query_labkey_rest_api, PROD_HOST

# Device and Assay Selection, to be incorporated into Lab Assist front end.
# Device: 'QuantStudio' uses a single .xls file (DupFile).
#         'LightCycler' uses one or two .txt files (DupFileFAM for RNase P, DupFileVIC for 16S).
# Assay:  '101' RNase P + 16S duplex, '61' RNase P only, '69' 16S only.

Device = 'QuantStudio'
Assay = '101'

# File Data, to be incorporated into Lab Assist front end.

ReagentFile = 'FRM-101-01_updated_draft.xlsx'
DupFile = 'SOP-101_v1_duplex_qPCR_day2_auto_QS1.xls'
DupFileFAM = 'day1_pipetmax_465.txt'
DupFileVIC = 'day1_pipetmax_533.txt'

# LabKey host used for every query, Lab Assist passes this in as host (defaults to production).
HOST = PROD_HOST

_BATCH_ID_RE = re.compile(r'^[A-Za-z0-9]+$')

def addsuffix(df: pd.DataFrame, col: str, group: str) -> pd.DataFrame:
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
        "/query-executeSql.api", host=host,
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
        ' CAST(ExtractedFrom.Name AS VARCHAR) AS ExtractedFrom'
        f' FROM samples."{assay_name}"'
        f" WHERE {batch_column} = '{batch_id}'"
        " ORDER BY Stored ASC"
    )
    result = query_labkey_rest_api(
        "/query-executeSql.api", host=host,
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
        "/query-executeSql.api", host=host,
        schemaName="core", sql=sql, columns="", includeHidden=True,
    )
    df = pd.DataFrame(result.get("rows", []))
    if df.empty:
        return pd.DataFrame(columns=['ExtractedFrom_1', 'Name'])
    return df[['ExtractedFrom_1', 'Name']]

# Import Reagent File to extract Run Information, use iloc because column names might change.

RunInfo = pd.read_excel(ReagentFile, header=3, usecols='B,E,F,G,I,J,K')
BatchID = RunInfo.iloc[0,4]
RunDate = RunInfo.iloc[0,5]
Operator = RunInfo.iloc[0,6]

# Asset ID depends on the thermocycler used, QuantStudio 1 and LightCycler 480 have their own rows.
if Device == 'QuantStudio':
    AssetID = str(RunInfo.iloc[14,2])
elif Device == 'LightCycler':
    AssetID = str(RunInfo.iloc[15,2])
else:
    raise ValueError(f"Device must be 'QuantStudio' or 'LightCycler', not '{Device}'")

# Import Reagent File again to get Lot Numbers, transpose for later use

Reagents = pd.read_excel(ReagentFile, header=3, usecols='B,E', nrows=8)
Reagents.set_index('Component', inplace=True)
Reagents.rename(columns={'Lot #':0}, inplace=True)
Reagents = Reagents.transpose()

if Device == 'QuantStudio':
    Results = pd.read_excel(DupFile, sheet_name='Results')


    # Find where Results actually start.
    FirstWell = Results.index[Results['Block Type'] == 'Well'][0]
    WellResults = pd.read_excel(DupFile, sheet_name='Results', header=FirstWell+1)
    WellResults['Sample Number'] = WellResults['Sample Name'].str.extract(r'(\d+)', expand=False).str.strip()
    WellResults['Sample Type'] = WellResults['Sample Name'].str.extract(r'(\D+)', expand=False).str.strip()

# Standard curve checks, shared by both devices.
def check_curve(SubAssay, Slope, RSqd, PCREff):

    ErrorString = []

    if pd.isna(Slope) or pd.isna(RSqd) or pd.isna(PCREff):
        ErrorString.append(f'\n{SubAssay} Standard Curve for {BatchID} is missing, concentrations and curve values will be blank.')
        return ErrorString

    if float(PCREff) < 85 or float(PCREff) > 115:
        ErrorString.append(f'\n{SubAssay} Standard Curve Efficiency for {BatchID} is outside of recommended range. Efficiency is recorded at {float(PCREff):.2f}%.')

    if float(RSqd) < 0.99:
        ErrorString.append(f'\n{SubAssay} Standard Curve R^2 for {BatchID} is outside of recommended range. R^2 is recorded at {float(RSqd):.4f}.')

    if float(Slope) > -3.10 or float(Slope) < -3.59:
        ErrorString.append(f'\n{SubAssay} Standard Curve Slope for {BatchID} is outside of recommended range. Standard Curve Slope is recorded at {float(Slope):.4f}.')

    return ErrorString

def MaleDNACheck(DNAVal):
    
    DNAString = []
    
    if not pd.isna(DNAVal):
        if DNAVal > 26.50 or DNAVal < 23.71:
            DNAString.append(f'\nHuman Male DNA for {BatchID} is outside of +/- 3 SD, with a value of {DNAVal:.3f}')
        elif DNAVal > 26.04 or DNAVal < 24.17:
            DNAString.append(f'\nHuman Male DNA for {BatchID} is outside of +/- 2 SD, with a value of {DNAVal:.3f}')
    return DNAString
       

# QuantStudio results, one .xls file holds both reporters.
def get_results(Reporter, SubAssay):
    GetStandards = WellResults.loc[~WellResults['Sample Type'].str.contains('Sample')].drop_duplicates(subset=['Sample Name','Reporter','Ct Mean', 'Ct SD']).copy()
    GetStandards = GetStandards[GetStandards['Reporter'] == Reporter].reset_index()
    Standards = GetStandards[['Sample Name', 'Ct Mean', 'Ct SD', ]].copy()
    Standards = Standards.set_index('Sample Name')
    # Human DNA Male is an RNase P control, so only check it on the FAM reporter.
    DNAString = []
    if Reporter == 'FAM':
        if 'Human DNA Male' in Standards.index:
            DNAMale = Standards.at['Human DNA Male','Ct Mean']
            DNAString = MaleDNACheck(DNAMale)
        else:
            DNAString = [f'\nHuman DNA Male control not found for {BatchID}, control check skipped.']
    Standards.rename(columns={'Ct Mean':Reporter + ' Mean Ct', 'Ct SD':Reporter + ' Ct SD'}, inplace=True)
    Standards = Standards.fillna('0')


    #Get curve values
    Slope = GetStandards['Slope'][0]
    YIntercept = GetStandards['Y-Intercept'][0]
    RSqd = GetStandards['R(superscript 2)'][0]
    PCREff = GetStandards['Efficiency'][0]

    ErrorString = check_curve(SubAssay, Slope, RSqd, PCREff)
    ErrorString = ''.join(ErrorString + DNAString)

    Standards2 = pd.DataFrame(data={0:[Slope, YIntercept, RSqd, PCREff]}, index=[SubAssay + ' Slope',SubAssay +  ' YIntercept',SubAssay +  ' RSqd',SubAssay +  ' PCREff'])
    Standards2 = Standards2.transpose()
    SampleWellResults = WellResults.loc[(WellResults['Sample Type'] != 'Standard')]
    SampleWellResults.loc[SampleWellResults['CT'] == 'Undetermined', 'CT'] = 0
    SampleWellResults = SampleWellResults[SampleWellResults['Reporter'] == Reporter]
    SampleWellResults = SampleWellResults[['Sample Name', 'Ct Mean', 'Ct SD', 'Sample Number', 'Sample Type', 'Quantity Mean']].drop_duplicates(subset='Sample Name')
    SampleWellResults.drop_duplicates(subset='Sample Name', inplace=True)
    SampleWellResults.dropna(subset='Sample Number', inplace=True)

    SampleWellResults = SampleWellResults[['Sample Name', 'Sample Number', 'Quantity Mean', 'Ct Mean', 'Ct SD']]
    SampleWellResults['Quantity Mean'] = SampleWellResults['Quantity Mean']/1000*20

    SampleWellResults.rename(columns={'Quantity Mean':SubAssay + ' Concentration (ng/μL)',
                           'Ct Mean':SubAssay + ' Mean CP',
                           'Ct SD':SubAssay + ' CP SD'},inplace=True)

    SampleWellResults['Sample Number'] = SampleWellResults['Sample Number'].astype(int)
    SampleWellResults.reset_index(inplace=True, drop=True)
    return SampleWellResults, Standards2, Standards, ErrorString

# LightCycler results, one .txt file per reporter.
def getLCResults(Filename, Reporter, SubAssay):    
    Results = pd.read_csv(Filename, sep='\t', header=1)

    # Well rows have a MeanCp, run statistics are the 'Name: value' lines after them.
    WellResults = Results[Results['MeanCp'].notna()].reset_index(drop=True)
    SampleOrder = pd.read_csv('Sample Order.csv')
    WellResults = pd.merge(WellResults, SampleOrder, how='inner', left_index=True, right_index=True)
    WellResults.rename(columns={'Sample ID':'Sample Name', 'MeanCp':'Ct Mean', 'STD Cp':'Ct SD'}, inplace=True)
    SampleWellResults = WellResults[WellResults['Sample Name'].str.contains('Sample')].copy()
    SampleWellResults['Sample Number'] = SampleWellResults['Sample Name'].str.extract(r'(\d+)', expand=False).str.strip().astype(int)
    SampleWellResults['Concentration (ng/uL)'] = SampleWellResults['Mean conc']/1000*20
    # LightCycler reports no amplification as Cp 0, blank these to match QuantStudio.
    SampleWellResults.loc[SampleWellResults['Ct Mean'] == 0, ['Ct Mean', 'Ct SD', 'Concentration (ng/uL)']] = math.nan
    SampleWellResults = SampleWellResults[['Sample Name', 'Sample Number','Concentration (ng/uL)','Ct Mean', 'Ct SD']]
    SampleWellResults.rename(columns={'Concentration (ng/uL)':SubAssay + ' Concentration (ng/μL)',
                                      'Ct Mean':SubAssay + ' Mean CP',
                                      'Ct SD':SubAssay + ' CP SD'},inplace=True)
    GetStandards = WellResults.loc[~WellResults['Sample Name'].str.contains('Sample')].copy()

    Standards = GetStandards[['Sample Name', 'Ct Mean', 'Ct SD', ]].copy()
    Standards = Standards.set_index('Sample Name')
    # Human DNA Male is an RNase P control, so only check it on the FAM reporter.
    DNAString = []
    if Reporter == 'FAM':
        if 'Human DNA Male' in Standards.index:
            DNAMale = Standards.at['Human DNA Male','Ct Mean']
            DNAString = MaleDNACheck(DNAMale)
        else:
            DNAString = [f'\nHuman DNA Male control not found for {BatchID}, control check skipped.']
    Standards.rename(columns={'Ct Mean':Reporter + ' Mean Ct', 'Ct SD':Reporter + ' Ct SD'}, inplace=True)
    
    Standards2 = Results[Results['Samples'].str.contains(': ', na=False)].copy()
    Standards2[['Type', 'Results']] = Standards2['Samples'].str.split(': ', expand=True)
    Standards2 = Standards2[['Type', 'Results']]
    Standards2['Results'] = Standards2['Results'].astype(float)
    Standards2.set_index('Type', inplace=True)
    Standards2 = Standards2.reindex(['Slope', 'YIntercept', 'Error', 'Efficiency'])
    Standards2.loc['Efficiency', 'Results'] = (10**(-1/Standards2.loc['Slope', 'Results']) - 1) * 100
    Standards2 = Standards2.transpose()
    Standards2['Error'] = 1-Standards2['Error']


    # Check for errors

    PCREff = Standards2.loc['Results', 'Efficiency']
    RSqd = Standards2.loc['Results', 'Error']
    Slope = Standards2.loc['Results', 'Slope']

    ErrorString = check_curve(SubAssay, Slope, RSqd, PCREff)
    ErrorString = ''.join(ErrorString + DNAString)

    Standards2.rename(columns={'Error':SubAssay + ' RSqd','Efficiency':SubAssay + ' PCREff', 'Slope':SubAssay + ' Slope', 'YIntercept':SubAssay + ' YIntercept'}, inplace=True)
    Standards2.columns.name = None
    SampleWellResults.reset_index(inplace=True, drop=True)
    Standards2.reset_index(inplace=True, drop=True)
    return SampleWellResults, Standards2, Standards, ErrorString

# Use the results function for the device selected, LightCycler has a .txt file per reporter.
def device_results(Reporter, SubAssay):
    if Device == 'QuantStudio':
        return get_results(Reporter, SubAssay)
    return getLCResults(DupFileFAM if Reporter == 'FAM' else DupFileVIC, Reporter, SubAssay)

if Assay == '101':

    RNaseP, RNasePStandards, RNasePSampleStandards, RNasePErrors = device_results('FAM', 'RNase P')
    RNaseP = pd.concat([RNaseP,RNasePStandards], axis=1)
    RNaseP[RNasePStandards.columns] = RNaseP[RNasePStandards.columns].ffill()

    PCR16S, PCR16SStandards, PCR16SSampleStandards, PCR16sErrors = device_results('VIC', '16S')
    PCR16S = pd.concat([PCR16S,PCR16SStandards], axis=1)
    PCR16S[PCR16SStandards.columns] = PCR16S[PCR16SStandards.columns].ffill()

    Standards = pd.concat([RNasePStandards, PCR16SStandards], axis=1)
    RunStandards = pd.merge(RNasePSampleStandards,PCR16SSampleStandards, how='inner', left_index=True, right_index=True)
    RunErrors = RNasePErrors + PCR16sErrors

    Results = pd.merge(RNaseP, PCR16S, on=['Sample Name','Sample Number'])
    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
           'RNase P Primer-Probe (FAM)',
           'Human Genomic DNA',
           'Human DNA Male control',
           '16S Primer-Probe (VIC)',
           'Nuclease free water or suitable alternative',
           'Microbial Community DNA Standard (2000ng)']]
    Results = pd.concat([Results, Reagents], axis=1)
    Results[Reagents.columns] = Results[Reagents.columns].ffill()
    Results['Asset ID'] = AssetID

    RunInfo = pd.DataFrame(data=[[Operator, BatchID, RunDate]], columns=['Operator', 'Batch ID', 'Processed Date'])

elif Assay == '61':

    Results, Standards, RunStandards, RunErrors = device_results('FAM', 'RNase P')
    Standards['Processed Date'] = RunDate

    RunInfo = pd.DataFrame(data=[[BatchID, Operator]], columns=['Batch ID', 'Operator'])

    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
           'RNase P Primer-Probe (FAM)',
           'Human Genomic DNA']]

    Reagents['Asset ID'] = AssetID

    RunInfo = pd.concat([RunInfo, Reagents, Standards], axis=1)
    RunInfo = RunInfo.transpose()

elif Assay == '69':

    Results, Standards, RunStandards, RunErrors = device_results('VIC', '16S')
    Standards['Processed Date'] = RunDate

    RunInfo = pd.DataFrame(data=[[BatchID, Operator]], columns=['Batch ID', 'Operator'])

    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
                         '16S Primer-Probe (VIC)',
           'Microbial Community DNA Standard (2000ng)']]

    Reagents['Asset ID'] = AssetID

    RunInfo = pd.concat([RunInfo, Reagents, Standards], axis=1)
    RunInfo = RunInfo.transpose()

# Batch ID is put into the LabKey SQL, so it must be alphanumeric.
batch_id = str(BatchID)
if not _BATCH_ID_RE.match(batch_id):
    raise ValueError(f"Invalid Batch ID '{batch_id}': must be alphanumeric (no spaces)")

# Fetch the manifest and the Human, Microbial and Bulk DNA for this batch from LabKey.
with ThreadPoolExecutor(max_workers=4) as executor:
    future_manifest = executor.submit(_fetch_manifest, batch_id, HOST)
    future_human = executor.submit(_fetch_dna_table, "Human DNA", "BatchID", batch_id, HOST)
    future_microbial = executor.submit(
        _fetch_dna_table, "Microbial DNA", "DEXID", batch_id, HOST,
    )
    future_bulk = executor.submit(_fetch_bulk_dna, batch_id, HOST)

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

dna = addsuffix(dna, 'ExtractedFrom_1', 'ExtractedFrom_1')
manifest = addsuffix(manifest, 'Sample_ID', 'Sample_ID')

manifest = manifest.merge(dna, how='inner', left_on='Sample_ID', right_on='ExtractedFrom_1')
manifest = manifest[['SampleNumber', 'Name']]

Results = Results.merge(manifest, how='inner', left_on='Sample Number', right_on='SampleNumber')

if Results.empty:
    raise ValueError(
        "Manifest and DNA were found, but no rows matched the Sample Numbers "
        "parsed from 'Sample Name' in the qPCR results."
    )

Results['Sample Name'] = Results['Name']
Results.drop(columns=['Sample Number','Name','SampleNumber'],inplace=True)

# Flag samples whose replicate CP SD is above 0.3, added to the same error log as the curve checks.
SubAssays = {'101': ['RNase P', '16S'], '61': ['RNase P'], '69': ['16S']}[Assay]
for SubAssay in SubAssays:
    HighSD = Results.loc[Results[SubAssay + ' CP SD'] > 0.3, ['Sample Name', SubAssay + ' CP SD']]
    for SampleName, SD in HighSD.itertuples(index=False):
        RunErrors += f'\n{SampleName} {SubAssay} CP SD for {BatchID} is above 0.3, recorded at {SD:.3f}.'

# Show Standard Curve, control and replicate errors and keep them with the output.
RunErrors = RunErrors.strip()
print(RunErrors if RunErrors else 'No run errors.')
RunErrors = pd.DataFrame({'Run Errors': RunErrors.split('\n') if RunErrors else ['None']})

with pd.ExcelWriter("Duplex Conversion.xlsx") as writer:
    Results.to_excel(writer, sheet_name="Results", index=False)
    Standards.to_excel(writer, sheet_name="Standards", index=False)
    RunStandards.to_excel(writer, sheet_name='Sample Standards')
    # 61 and 69 Run Information is transposed, so its labels are in the index.
    RunInfo.to_excel(writer, sheet_name='Run Information', index=(Assay != '101'))
    RunErrors.to_excel(writer, sheet_name='Run Errors', index=False)

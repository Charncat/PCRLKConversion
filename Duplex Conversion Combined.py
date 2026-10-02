# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 2026

@author: DanielMoutia

Combined Duplex Conversion for QuantStudio 1 (.xls) and LightCycler 480 (.txt) exports.
Merges Duplex Conversion Excel.py and Duplex Conversion Notepad.py.
"""

# Import all packages (remove any that are not needed)

import pandas as pd
import math
import requests
import os
import sys
import json
from concurrent.futures import ThreadPoolExecutor

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

BASE   = "https://originsciences.app.labkey.host"
FOLDER = "home"
SCHEMA = "samples"
# API key is read from the environment so it is not stored in the code.
APIKEY = os.environ["LABKEY_API_KEY"]

s = requests.Session()
s.headers["Authorization"] = f"LABKEY apikey={APIKEY}"

def query_labkey(path, **params):
    r = s.get(f"{BASE}/home/{path}", params={"containerPath": FOLDER, "columns": "", **params}, timeout=30)

    try:
        r.raise_for_status()
    except requests.HTTPError:
        print("Status:", r.status_code)
        try:
            print(r.json(),)
        except Exception:
            print(r.text)

    r.raise_for_status()
    return r.json()

def query_labkey_rest_api(path_controller_action: str, **params: str) -> json:
  host = "originsciences.app.labkey.host"
  container = "home"
  try:
    request = requests.get(
      f"https://{host}/home/{path_controller_action}",
      params={"containerPath": container, **params}
    )
  except requests.HTTPError as e:
    print(f"Error: HTTP request to LabKey failed: {path_controller_action}")
    print("This could be due to invalid credentials or a lack of permissions")
    print(f"{type(e).__name__} was raised: {e}")
    sys.exit(1)
  return request.json()

def find_templates() -> list:
  queries = query_labkey_rest_api("/query-GetQueries.api", schemaName="workflow.jobtypes", maxRows=None).get("queries", [])
  templates_list = [x.get("name") for x in queries]
  return templates_list

def query_workflow_fields(batch_id: str) -> pd.DataFrame:
  job_templates = find_templates()

  def fetch_template(job_template):
    sql = f"""
    SELECT
    jt.Name,
    jt.BatchID,
    CAST(jt.QubitOperator.DisplayName AS VARCHAR) AS QubitOperator
    FROM {job_template} as jt
    WHERE jt.BatchID = '{batch_id}'
    """
    return pd.DataFrame(query_labkey_rest_api("/query-executeSql.api", schemaName="workflow.jobtypes",
                                              sql=sql, maxRows=None).get("rows", []))

  with ThreadPoolExecutor() as executor:
    dfs = list(executor.map(fetch_template, job_templates))

  workflows_df = pd.concat(dfs, axis=0, ignore_index=True) if dfs else pd.DataFrame(columns=["Name", "BatchID", "QubitOperator"])
  return workflows_df[["Name", "BatchID", "QubitOperator"]]

def getjobs():
  df = query_workflow_fields(BatchID)
  df.fillna('-', inplace=True)
  #df = df[~df.isnull().any(axis=1)]
  return df

# Function to deduplicate rows by adding a suffix
def addsuffix(df, col, group):

    df[col] += df.groupby(group).cumcount().add(1).astype(str).radd('_').mask(df.groupby(group)[col].transform('count')==1,'')
    return df

# Function to pull all Human DNA and what OriCol Aliquots they were extracted from.
def get_samples(assay_name, BatchID):
    sql=f'SELECT Name,ExtractedFrom,{BatchID} AS BatchID, CAST(ExtractedFrom.Name AS VARCHAR) AS ExtractedFrom FROM samples."{assay_name}"'

    results = pd.DataFrame(query_labkey("/query-executeSql.api", schemaName="core", sql=sql, maxRows=max_rows, includeHidden=True).get("rows", []))
    return results

def get_samples2(assay_name):
    sql=f'SELECT Name,SecondaryParent,ExperimentID AS BatchID, CAST(SecondaryParent.Name AS VARCHAR) AS ExtractedFrom_1 FROM samples."{assay_name}"'

    results = pd.DataFrame(query_labkey("/query-executeSql.api", schemaName="core", sql=sql, maxRows=max_rows, includeHidden=True).get("rows", []))
    return results

# Fire all three API fetches concurrently — they only depend on BatchID (already extracted above).
custom_sql_select_query = ", CAST(SampleID.Name AS VARCHAR) AS Sample_ID, CAST(Run.BatchID AS VARCHAR) AS Batch_ID"

def _fetch_manifest(assay_name):
    manifest_sql = f'SELECT *{custom_sql_select_query} FROM assay.General."{assay_name}".Data AS d'
    return pd.DataFrame(query_labkey("/query-executeSql.api", schemaName="core", sql=manifest_sql,
                                     maxRows=max_rows, includeHidden=True).get("rows", []))

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

max_rows = None
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

with ThreadPoolExecutor(max_workers=5) as executor:
    future_manifest_98 = executor.submit(_fetch_manifest, "SOP-98 Manifest Creation")
    future_manifest_64 = executor.submit(_fetch_manifest, "SOP-64 Manifest Creation")
    future_human      = executor.submit(get_samples, "Human DNA", "BatchID")
    future_bulk     = executor.submit(get_samples2, "DNA")
    future_microbial = executor.submit(get_samples, "Microbial DNA", "DEXID")

# API call to get manifest, contains list of samples and their 'numbers', filtered down to relevant Run ID and then positions converted to just number.
Manifest = pd.concat([future_manifest_98.result(), future_manifest_64.result()], ignore_index=True)
Manifest = Manifest[['ProcessPosition', 'Batch_ID', 'Sample_ID']]
Manifest = Manifest[Manifest['Batch_ID'] == BatchID]
Manifest['SampleNumber'] = Manifest['ProcessPosition'].str.extract(r'(\d+)').astype(int)
Manifest = Manifest[['Sample_ID', 'SampleNumber']]

# All Human DNA is called and then filtered down to the relevant Run ID
DNAHu = future_human.result()
DNAMb = future_microbial.result()
DNABk = future_bulk.result()
DNAHu = pd.concat([DNAHu, DNAMb, DNABk])
DNAHu = DNAHu[DNAHu['BatchID'] == BatchID]
DNAHu = DNAHu[['ExtractedFrom_1', 'Name']]


# Dedupe both Human DNA and Manifest.
DNAHu = addsuffix(DNAHu,'ExtractedFrom_1','ExtractedFrom_1')
Manifest = addsuffix(Manifest,'Sample_ID','Sample_ID')

# Merge both dataframes, drop uneeded columns.
Manifest = Manifest.merge(DNAHu, how='inner', left_on='Sample_ID', right_on='ExtractedFrom_1')
Manifest = Manifest[['SampleNumber', 'Name']]

Results = Results.merge(Manifest, how='inner', left_on='Sample Number', right_on='SampleNumber')
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

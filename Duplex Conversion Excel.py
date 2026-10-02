# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:31:09 2026

@author: DanielMoutia
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:32:47 2026

@author: DanielMoutia
"""

import re
import pandas as pd
import io
import math
import requests
import os
import sys
import json
from concurrent.futures import ThreadPoolExecutor


DupFile = 'SOP-101_v1_duplex_qPCR_day2_manual_QS1.xls'
ReagentFile = 'FRM-101-01_updated_draft.xlsx'
Assay = '101'

BASE   = "https://originsciences.app.labkey.host"
FOLDER = "home"
SCHEMA = "samples"
APIKEY = "5830a18d79e74c1589bda8df43762f1a348d0ad85d243c075465df6329fa4b29"

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

# Name search function
def namesearch(regex):
    for filename in os.listdir():
        if re.search(regex, filename):
            return filename

# Import Reagent File to extract Run Information, use iloc because column names might change.

RunInfo = pd.read_excel(ReagentFile, header=3, usecols='B,E,F,G,I,J,K') 
PConc = RunInfo.loc[6, 'Conc, ng/uL']
PUsed = RunInfo.loc[6, 'uL for 60uL STD#1']
BatchID = RunInfo.iloc[0,4]
RunDate = RunInfo.iloc[0,5]
Operator = RunInfo.iloc[0,6]
AssetID = str(RunInfo.iloc[15,2])

# Import Reagent File again to get Lot Numbers, transpose for later use

Reagents = pd.read_excel(ReagentFile, header=3, usecols='B,E', nrows=8) 
Reagents.set_index('Component', inplace=True)
Reagents.rename(columns={'Lot #':0}, inplace=True)
Reagents = Reagents.transpose()

max_rows = None
if DupFile[-3:] == 'xls':
    Results = pd.read_excel(DupFile, sheet_name='Results')

  
    # Find where Results actually start.
    FirstWell = Results.index[Results['Block Type'] == 'Well'][0]
    WellResults = pd.read_excel(DupFile, sheet_name='Results', header=FirstWell+1)
    WellResults['Sample Number'] = WellResults['Sample Name'].str.extract(r'(\d+)', expand=False).str.strip()
    WellResults['Sample Type'] = WellResults['Sample Name'].str.extract(r'(\D+)', expand=False).str.strip()
   
def get_results(Reporter, SubAssay):
    GetStandards = WellResults.loc[~WellResults['Sample Type'].str.contains('Sample')].drop_duplicates(subset=['Sample Name','Reporter','Ct Mean', 'Ct SD']).copy()
    GetStandards = GetStandards[GetStandards['Reporter'] == Reporter].reset_index()
    Standards = GetStandards[['Sample Name', 'Ct Mean', 'Ct SD', ]].copy()
    Standards['Sample Name'] = SubAssay + ' ' + Standards['Sample Name']
    Standards = Standards.set_index('Sample Name').stack(future_stack=True).to_frame().T
    Standards.columns = [f"{name} {metric}" for name, metric in Standards.columns]
    Standards = Standards.transpose()
    Standards = Standards.fillna('0')

    
    #Get curve values    
    Slope = GetStandards['Slope'][0]
    YIntercept = GetStandards['Y-Intercept'][0]
    RSqd = GetStandards['R(superscript 2)'][0]
    PCREff = GetStandards['Efficiency'][0]

    ErrorString = []
    
    if float(PCREff) < 85 or float(PCREff) > 115:
        ErrorString.append(f'\nStandard Curve Efficiency for {BatchID} is outside of recommended range. Effiency is recorded at {PCREff}%.\n')
        
    if RSqd < 0.99:
        ErrorString.append(f'\n Standard Curve R^2 for {BatchID} is outside of recommended range. R^2 is recorded at {RSqd}')
        
    if float(Slope) > -3.10 or float(Slope) < -3.59:
        ErrorString.append(f'\n Standard Curve Slope for {BatchID} is outside of recommdend range. Standard Curve Slope is recorded at {Slope}')
 
    ErrorString = ''.join(ErrorString)
    
    Standards2 = pd.DataFrame(data={0:[Slope, YIntercept, RSqd, PCREff]}, index=[SubAssay + ' Slope',SubAssay +  ' YIntercept',SubAssay +  ' RSqd',SubAssay +  ' PCREff'])
    Standards2 = Standards2.transpose()
    SampleWellResults = WellResults.loc[(WellResults['Sample Type'] != 'Standard')]
    SampleWellResults.loc[SampleWellResults['CT'] == 'Undetermined', 'CT'] = 0
    SampleWellResults = SampleWellResults[SampleWellResults['Reporter'] == Reporter]   
    SampleWellResults = SampleWellResults[['Sample Name', 'Ct Mean', 'Ct SD', 'Sample Number', 'Sample Type', 'Quantity Mean']].drop_duplicates(subset='Sample Name')
    SampleWellResults.drop_duplicates(subset='Sample Name', inplace=True)
    SampleWellResults.dropna(subset='Sample Number', inplace=True)
    
    SampleWellResults = SampleWellResults[['Sample Name', 'Sample Number', 'Quantity Mean', 'Ct Mean', 'Ct SD']]
    
    SampleWellResults.rename(columns={'Quantity Mean':SubAssay + ' Concentration (ng/μL)', 
                           'Ct Mean':SubAssay + ' Mean CP ', 
                           'Ct SD':SubAssay + ' CP SEM'},inplace=True)
    
    SampleWellResults['Sample Number'] = SampleWellResults['Sample Number'].astype(int)
    SampleWellResults.reset_index(inplace=True, drop=True)
    return SampleWellResults, Standards2, Standards, ErrorString

if Assay == '101':

    RNaseP, Standards, RNasePSampleStandards, RNasePErrors = get_results('FAM', 'RNase P')
    RNaseP = pd.concat([RNaseP,Standards], axis=1)
    RNaseP.ffill(inplace=True)
    
    PCR16S, Standards, PCR16SSampleStandards, PCR16sErrors = get_results('VIC', '16s')
    PCR16S = pd.concat([PCR16S,Standards], axis=1)
    PCR16S.ffill(inplace=True)
    
    RunStandards = pd.concat([RNasePSampleStandards,PCR16SSampleStandards])
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
    Results.ffill(inplace=True)
    Results['Asset ID'] = AssetID     

    RunInfo = pd.DataFrame(data=[[Operator, BatchID, RunDate]], columns=['Operator', 'Batch ID', 'Processed Date'])

elif Assay == '61':
    
    Results, Standards, RunStandards, RunErrors = get_results('FAM', 'RNase P')
    Standards = Standards.transpose()
    Standards['Processed Date'] = RunDate
    
    RunInfo = pd.DataFrame(data=[[BatchID, Operator]], columns=['Batch ID', 'Operator'])
        
    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
           'RNase P Primer-Probe (FAM)', 
           'Human Genomic DNA']]
    
    Reagents['Asset ID'] = AssetID
    
    RunInfo = pd.concat([RunInfo, Reagents, Standards], axis=1)
    RunInfo = RunInfo.transpose()

elif Assay == '69':

    Results, Standards, RunStandards, RunErrors = get_results('VIC', '16S')
    Standards = Standards.transpose()    
    Standards['Processed Date'] = RunDate
    
    RunInfo = pd.DataFrame(data=[[BatchID, Operator]], columns=['Batch ID', 'Operator']) 
    
    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
                         '16S Primer-Probe (VIC)',
           'Microbial Community DNA Standard (2000ng)']]
    
    Reagents['Asset ID'] = AssetID
    
    RunInfo = pd.concat([RunInfo, Reagents, Standards], axis=1)
    RunInfo = RunInfo.transpose()

with ThreadPoolExecutor(max_workers=6) as executor:
    future_manifest_98 = executor.submit(_fetch_manifest, "SOP-98 Manifest Creation")
    future_manifest_64 = executor.submit(_fetch_manifest, "SOP-64 Manifest Creation")
    future_human      = executor.submit(get_samples, "Human DNA", "BatchID")
    future_bulk     = executor.submit(get_samples2, "DNA")
    future_microbial = executor.submit(get_samples, "Microbial DNA", "DEXID")
    future_jobs     = executor.submit(getjobs)

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
Results.drop(columns=['SampleNumber','Name'],inplace=True)
    
# with pd.ExcelWriter("Duplex Conversion.xlsx") as writer:
#     Results.to_excel(writer, sheet_name="Results", index=False)
#     Standards.to_excel(writer, sheet_name="Standards", index=False)
#     RunStandards.to_excel(writer, sheet_name='Sample Standards')
#     RunInfo.to_excel(writer, sheet_name='Run Information', index=False)
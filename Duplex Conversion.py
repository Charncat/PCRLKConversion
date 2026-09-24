# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 09:01:59 2026

@author: DanielMoutia
"""

DupFile = '20251105_TRIOMIC_RNaseP-16S_Hu_DEN078.txt'

### EXTERNAL VARIABLES TO BE ENTERED BY USER ###

Promega = 209
ExpectedMeanConcs = [20000, 2500, 312.50, 39.06, 4.88, 0.61]

### CALCULATED VARIABLES

PromegaUsed = round(150000 / (Promega*1000) * 20, 2)
if int(str(PromegaUsed)[-2:]) % 2 == 1:
    PromegaUsed += 0.01

import re
import pandas as pd
import io
import math
import requests
import os
import sys
import json
from concurrent.futures import ThreadPoolExecutor

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
  queries = query_labkey_rest_api("/query-GetQueries.api", schemaName="sampleManagement.jobtypes", maxRows=None).get("queries", [])
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
    return pd.DataFrame(query_labkey_rest_api("/query-executeSql.api", schemaName="sampleManagement.jobtypes",
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

# Name search function
def namesearch(regex):
    for filename in os.listdir():
        if re.search(regex, filename):
            return filename
max_rows = None


text = open(DupFile).read()
#Run = re.search(r'Experiment Name = \d+-\d+-\d+_\d+_(\D+\d+)_\D+-\d+\D', text)[1]
BatchID = 'DEN078'
result = re.findall(r'\[Results\](.+)', text, flags=re.DOTALL)[0]
PCR = pd.read_table(io.StringIO(result))
PCR['Sample Number'] = PCR['Sample Name'].str.extract(r'(\d+)', expand=False).str.strip()
PCR['Sample Type'] = PCR['Sample Name'].str.extract(r'(\D+)', expand=False).str.strip()

#Calculate Mean CPs, Mean Concentration, and Log Concentration

#Get Means and SEMs of RNase Standards
def get_results(Reporter, SubAssay):
    Slope = PCR.loc[~PCR['Sample Type'].str.contains('Sample')].drop_duplicates(subset='Ct Mean').copy()
    Slope = Slope[Slope['Reporter'] == Reporter]
    Standards = Slope[['Sample Name', 'Ct Mean', 'Ct SD']].copy()
    Standards['Sample Name'] = SubAssay + ' ' + Standards['Sample Name']
    Standards = Standards.set_index('Sample Name').stack(dropna=False).to_frame().T
    Standards.columns = [f"{name} {metric}" for name, metric in Standards.columns]
    Standards = Standards.transpose()
    
    Slope = PCR.loc[(PCR['Sample Type'] == 'Standard')].drop_duplicates(subset='Ct Mean')
    Slope = Slope[Slope['Reporter'] == Reporter]
    Slope = Slope[['Sample Name', 'Ct Mean']]
    AMeanCp = Slope['Ct Mean'].tolist()
    
    EMeanCp = []
    MeanConc = Promega*1000 * PromegaUsed/20
    MeanConcList = []
    for i in AMeanCp:
        MeanConcList.append(MeanConc)
        EMeanCp.append(math.log(((3.6*(10**11))/((2*(MeanConc/1000))/0.00646246)),2))
        MeanConc = MeanConc*(4.5/20)
      
    LogConc = []
    for i in MeanConcList:
        LogConc.append(math.log(i, 10))
    
    #Get curve values    
    from scipy.stats import linregress
    Slope = linregress(LogConc,AMeanCp)[0]
    YIntercept = linregress(LogConc,AMeanCp)[1]
    RSqd = (linregress(LogConc,AMeanCp)[2])**2
    PCREff = ((10**(-1/Slope))-1)*100
    
    Standards2 = pd.DataFrame(data={0:[Slope, YIntercept, RSqd, PCREff]}, index=[SubAssay + ' Slope',SubAssay +  ' YIntercept',SubAssay +  ' RSqd',SubAssay +  ' PCREff'])
    AllStandards = pd.concat([Standards,Standards2])
    # Get mean Concentration
    AMeanConc = []
    for i in AMeanCp:
        AMeanConc.append(10**((i - YIntercept)/Slope))
        
    PCRExp = PCR.loc[(PCR['Sample Type'] != 'Standard')]
    PCRExp.loc[PCR['CT'] == 'Undetermined', 'CT'] = 0
    PCRExp = PCRExp[PCRExp['Reporter'] == 'FAM']
    PCRExpSEM = PCRExp.groupby(['Sample Name'])
    GroupedA = PCRExpSEM['CT'].std(ddof=0)
    GroupedA = GroupedA/math.sqrt(3)
    
    PCRExp = PCRExp[['Sample Name', 'Ct Mean', 'Ct SD', 'Sample Number', 'Sample Type']].drop_duplicates(subset='Sample Name')
    PCRExp[['Slope', 'YIntercept', 'Rsquared', 'PCREfficiency']] = [Slope, YIntercept, RSqd, PCREff]
    PCRExp['SampleConcentration'] = 10**((PCRExp['Ct Mean'] - YIntercept)/Slope)
    PCRExp['SampleConcentration'] = PCRExp['SampleConcentration']/100
    PCRExp.drop_duplicates(subset='Sample Name', inplace=True)
    PCRExp = PCRExp.merge(GroupedA, left_on='Sample Name', right_index=True)
    PCRExp.dropna(subset='Sample Number', inplace=True)
    
    PCRExp = PCRExp[['Sample Name', 'Sample Number', 'SampleConcentration', 'Ct Mean', 'CT']]
    
    PCRExp.rename(columns={'SampleConcentration':SubAssay + ' Concentration (ng/μL)', 
                           'Ct Mean':SubAssay + ' Mean CP ', 
                           'CT':SubAssay + ' CP SEM'},inplace=True)
    
    PCRExp['Sample Number'] = PCRExp['Sample Number'].astype(int)
    
    return PCRExp, AllStandards

PCRRNase, StandardsRNase = get_results('FAM', 'RNase P')
PCR16S, Standards16S = get_results('VIC', '16S')

# LogConc = []
# for i in ExpectedMeanConcs:
#     LogConc.append(math.log(i, 10))
    
# Slope = PCR.loc[(PCR['Sample Type'] == 'Standard')].drop_duplicates(subset='Ct Mean')
# Slope = Slope[Slope['Reporter'] == 'VIC']
# Slope = Slope[['Sample Name', 'Ct Mean']]
# AMeanCp = Slope['Ct Mean'].tolist()


# #Get curve values    
# from scipy.stats import linregress
# Slope = linregress(LogConc,AMeanCp)[0]
# YIntercept = linregress(LogConc,AMeanCp)[1]
# RSqd = (linregress(LogConc,AMeanCp)[2])**2
# PCREff = ((10**(-1/Slope))-1)*100

# # Get mean Concentration
# AMeanConc = []
# for i in AMeanCp:
#     AMeanConc.append(10**((i - YIntercept)/Slope))
    
# PCRExp = PCR.loc[(PCR['Sample Type'] != 'Standard')]
# PCRExp.loc[PCR['CT'] == 'Undetermined', 'CT'] = 0
# PCRExp = PCRExp[PCRExp['Reporter'] == 'VIC']
# PCRExpSEM = PCRExp.groupby(['Sample Name'])
# GroupedA = PCRExpSEM['CT'].std(ddof=0)
# GroupedA = GroupedA/math.sqrt(3)

# PCRExp = PCRExp[['Sample Name', 'Ct Mean', 'Ct SD', 'Sample Number', 'Sample Type']].drop_duplicates(subset='Sample Name')
# PCRExp[['Slope', 'YIntercept', 'Rsquared', 'PCREfficiency']] = [Slope, YIntercept, RSqd, PCREff]
# PCRExp['SampleRNaseConcentration'] = 10**((PCRExp['Ct Mean'] - YIntercept)/Slope)
# PCRExp.drop_duplicates(subset='Sample Name', inplace=True)
# PCRExp = PCRExp.merge(GroupedA, left_on='Sample Name', right_index=True)
# PCRExp['SampleRNaseConcentration'] = PCRExp['SampleRNaseConcentration']/100
# PCRExp.dropna(subset='Sample Number', inplace=True)

                    
# PCRExp = PCRExp[['Sample Name', 'Sample Number', 'SampleRNaseConcentration', 'Ct Mean', 'CT',
#         'Rsquared', 'PCREfficiency', 'YIntercept', 'Slope']]

# PCRExp.rename(columns={'SampleRNaseConcentration':'16S Concentration (ng/μL)', 
#                        'Ct Mean':'16S Mean CP ', 
#                        'CT':'16S CP SEM',
#                     'Rsquared':'16S R^2', 
#                     'PCREfficiency':'16S % PCR Efficiency', 
#                     'YIntercept':'16S Y-Intercept', 
#                     'Slope':'16S Slope'},inplace=True)

# PCRExp['Sample Number'] = PCRExp['Sample Number'].astype(int)
Standards = pd.concat([Standards16S, StandardsRNase])
Duplex = pd.merge(PCRRNase, PCR16S, how='left', on=['Sample Number', 'Sample Name'])
Standards.to_csv('Standards.csv')
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

Duplex = Duplex.merge(Manifest, how='inner', left_on='Sample Number', right_on='SampleNumber')
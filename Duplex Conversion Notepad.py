# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:50:22 2026

@author: DanielMoutia
"""

# Import all packages (remove any that are not needed)

import re
import pandas as pd
import io
import math
import requests
import os
import sys
import json
from concurrent.futures import ThreadPoolExecutor
from scipy.stats import linregress

# Assay Selection, to be incorporated into front end.

Assay = '101'

# File Data, to be incorporated into front end.

ReagentFile = 'FRM-101-01_updated_draft.xlsx'
DupFileFAM = 'LC480_pipetmax_465_derivative_max.txt'
DupFileVIC = 'LC480_pipetmax_533_derivative_max.txt'

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

def getLCResults(Filename, SubAssay):
    Filename = 'LC480_manual_465_derivative_max.txt'
    SubAssay = 'RNase P'
    Results = pd.read_csv(Filename, sep='\t', header=1)

    WellResults = Results[:-5]
    SampleOrder = pd.read_csv('Sample Order.csv')
    WellResults = pd.merge(WellResults, SampleOrder, how='inner', left_index=True, right_index=True)
    SampleWellResults = WellResults[WellResults['Sample ID'].str.contains('Sample')]
    SampleWellResults['Sample Number'] = SampleWellResults['Sample ID'].str.extract(r'(\d+)', expand=False).str.strip()
    SampleWellResults['Concentration (ng/uL)'] = SampleWellResults['Mean conc']/100
    SampleWellResults = SampleWellResults[['Sample ID', 'Sample Number','Concentration (ng/uL)','MeanCp', 'STD Cp']]
    SampleWellResults.rename(columns={'Concentration (ng/uL)':SubAssay + 'Concentration (ng/uL)', 
                                      'MeanCp':SubAssay + ' Mean CP',
                                      'STD Cp':SubAssay + ' CP SEM'},inplace=True)
    GetStandards = WellResults.loc[~WellResults['Sample ID'].str.contains('Sample')].copy()
    
    Standards = GetStandards[['Sample ID', 'MeanCp', 'STD Cp', ]].copy()
    Standards['Sample Name'] = SubAssay + ' ' + Standards['Sample ID']
    Standards = Standards.set_index('Sample Name').stack(future_stack=True).to_frame().T
    Standards.columns = [f"{name} {metric}" for name, metric in Standards.columns]
    Standards = Standards.transpose()
    Standards = Standards.fillna('0')
    
    Standards2 = Results[-5:].copy()
    Standards2[['Type', 'Results']] = Standards2['Samples'].str.split(': ', expand=True)
    Standards2 = Standards2[['Type', 'Results']]
    Standards2.iloc[1,1] = (float(Standards2.iloc[1,1])-1)*100
    Standards2.set_index('Type', inplace=True)
    Standards2 = Standards2.transpose()
    Standards2['Error'] = 1-float(Standards2['Error'])
    Standards2 = Standards2[['Error', 'Efficiency', 'YIntercept', 'Slope']]
    
    
    # Check for errors
    
    PCREff = float(Standards2.loc['Results', 'Efficiency'])
    RSqd = float(Standards2.loc['Results', 'Error'])
    Slope = float(Standards2.loc['Results', 'Slope'])
    
    ErrorString = []
    
    if PCREff < 85 or PCREff > 115:
        ErrorString.append(f'\nStandard Curve Efficiency for {BatchID} is outside of recommended range. Effiency is recorded at {PCREff}%.\n')
        
    if RSqd < 0.99:
        ErrorString.append(f'\n Standard Curve R^2 for {BatchID} is outside of recommended range. R^2 is recorded at {RSqd}')
        
    if Slope > -3.10 or Slope < -3.59:
        ErrorString.append(f'\n Standard Curve Slope for {BatchID} is outside of recommdend range. Standard Curve Slope is recorded at {Slope}')
    
    ErrorString = ''.join(ErrorString)

    Standards2.rename(columns={'Error':SubAssay + ' R^2','Efficiency':SubAssay + ' % PCR Efficiency', 'Slope':SubAssay + ' Slope', 'YIntercept':SubAssay + ' Y-Intercept'}, inplace=True)
    SampleWellResults.reset_index(inplace=True, drop=True)
    Standards2.reset_index(inplace=True, drop=True)
    return SampleWellResults, Standards2, ErrorString, Standards

if Assay == '101':

    RNaseP, Standards, RNasePErrors, RNasePRunStandards = getLCResults(DupFileFAM, 'RNase P')
    RNaseP = pd.concat([RNaseP,Standards], axis=1)
    RNaseP.ffill(inplace=True)
    
    PCR16S, Standards, PCR16sErrors, PCR16sRunStandards = getLCResults(DupFileVIC, '16s')
    PCR16S = pd.concat([PCR16S,Standards], axis=1)
    PCR16S.ffill(inplace=True)
    
    RunErrors = RNasePErrors + PCR16sErrors
    RunStandards = pd.concat([RNasePRunStandards, PCR16sRunStandards])
    
    Results = pd.merge(RNaseP, PCR16S, on=['Sample ID','Sample Number'])
    Results = Results.apply(pd.to_numeric, errors='ignore')
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

    Results, Standards, ErrorString, RunStandards = getLCResults(DupFileFAM, 'RNase P')
    Standards = Standards.apply(pd.to_numeric, errors='ignore')
    Standards['Processed Date'] = RunDate
    
    RunInfo = pd.DataFrame(data=[[BatchID, Operator]], columns=['Batch ID', 'Operator'])
        
    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
           'RNase P Primer-Probe (FAM)', 
           'Human Genomic DNA']]
    
    Reagents['Asset ID'] = AssetID
    
    RunInfo = pd.concat([RunInfo, Reagents, Standards], axis=1)
    RunInfo = RunInfo.transpose()
    
elif Assay == '69':

    Results, Standards, ErrorString, RunStandards = getLCResults(DupFileVIC, '16s')
    Standards = Standards.apply(pd.to_numeric, errors='ignore')
    Standards['Processed Date'] = RunDate
    
    RunInfo = pd.DataFrame(data=[[BatchID, Operator]], columns=['Batch ID', 'Operator']) 
    
    Reagents = Reagents[['TaqMan® Fast Advanced Master Mix (5x5mL)',
                         '16S Primer-Probe (VIC)',
           'Microbial Community DNA Standard (2000ng)']]
    
    Reagents['Asset ID'] = AssetID
    
    RunInfo = pd.concat([RunInfo, Reagents, Standards], axis=1)
    RunInfo = RunInfo.transpose()
    
with pd.ExcelWriter("Duplex Conversion.xlsx") as writer:
    Results.to_excel(writer, sheet_name="Results", index=False)
    Standards.to_excel(writer, sheet_name="Standards", index=False)
    RunStandards.to_excel(writer, sheet_name='Sample Standards')
    RunInfo.to_excel(writer, sheet_name='Run Information', index=False)

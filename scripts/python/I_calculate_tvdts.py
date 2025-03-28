"""
Script to calculate tumour volume doubling times for the 34 cases
"""
import pandas as pd
import numpy as np
import sys
import os


# Paths
path_to_volumes = './../../data/raw_volumes.csv'
path_for_output = './../../output/'
os.makedirs(path_for_output, exist_ok=True)

vols = pd.read_csv(path_to_volumes)

ratios_ov = vols[vols.valid_ov.ge(1)].groupby('anon_id').vol_ov.apply(lambda x: x.iloc[-1]/ x.iloc[0])
ratios_ov.name = 'ratio_ov'
ratios_om = vols[vols.valid_om.ge(1)].groupby('anon_id').vol_om.apply(lambda x: x.iloc[-1]/ x.iloc[0])
ratios_om.name = 'ratio_om'
dt = vols.groupby('anon_id').dt.apply(lambda x: x.iloc[-1])
dt.name = 'dt'

'''
V = V0 exp(r*t)
r = 1/t ln(V/V0)
tvdt = 1/r ln(2) = t * ln(2) / ln(V/V0)
'''
tvdt = pd.concat([dt, ratios_ov, ratios_om], axis=1).reset_index()
tvdt['tvdt_ov'] = tvdt.apply(lambda x: x.loc['dt'] * np.log(2) / np.log(x.loc['ratio_ov']), axis=1)
tvdt['tvdt_om'] = tvdt.apply(lambda x: x.loc['dt'] * np.log(2) / np.log(x.loc['ratio_om']), axis=1)

tvdt.to_csv(path_for_output + 'tvdts.csv')
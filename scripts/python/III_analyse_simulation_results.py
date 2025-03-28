'''
Script to take data from simulations and determine
1. What proportion of tumours will metastasise before the earliest possible screen detection
2. What the window of opportunity is for those tumours that can be detected before metastasis.
'''
import pandas as pd
import os
import matplotlib.pyplot as plt


# Paths
path_to_sims = './../../output/simulations/sims.csv'

sims = pd.read_csv(path_to_sims, index_col=0)

'''
Statistics
1a. How many reach CA125 limit before mets
1b. What is the WOO between CA125 limit and mets
2a. How many reach US detection limit before mets
2b. What is the WOO between US limit and mets
'''

ca125_b4_mets = sims[sims.time_to_met.gt(sims.time_to_ca125)]
US_b4_mets = sims[sims.time_to_met.gt(sims.time_to_US)]


print('{} out of 10,000 cases reach CA125 detection limit before mets'.format(sims[sims.time_to_met.gt(sims.time_to_ca125)].shape[0]))
print('{} out of 10,000 cases reach US detection limit before mets'.format(sims[sims.time_to_met.gt(sims.time_to_US)].shape[0]))

print('\n **** WOO for CA125 stats **** \n')
print((ca125_b4_mets.time_to_met * 12 / 365 - ca125_b4_mets.time_to_ca125 * 12 / 365).agg('describe'))
print('\n')

print('\n **** WOO for US stats **** \n')
print((US_b4_mets.time_to_met * 12 / 365 - US_b4_mets.time_to_US * 12 / 365).agg('describe'))
print('\n')

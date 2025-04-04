'''
TVDT plotting for the 34 cases 
'''
import pandas
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


# Paths
path_to_volumes = '../../../data/raw_volumes.csv'
path_to_tvdts = '../../../output/tvdts.csv'
path_for_output = '../../../output/plots/'

if not os.path.exists(path_for_output):
    os.makedirs(path_for_output)

volumes = pd.read_csv(path_to_volumes)
volumes.sort_values(by=['anon_id', 'dt'])
growth_rates = pd.read_csv(path_to_tvdts, index_col=0)

'''
Statistics
'''
print('\n **** Stats TVDT all **** \n')
print('## Ovarian \n')
print(growth_rates[growth_rates.tvdt_ov.gt(0)].tvdt_ov.apply(lambda x: x * 12 / 365).agg('describe')) # In months
print('## \n Omental \n')
print(growth_rates[growth_rates.tvdt_om.gt(0)].tvdt_om.apply(lambda x: x * 12 / 365).agg('describe')) # In months


'''
1. Box plots for all omental and ovarian lesions together
'''
X = [growth_rates[growth_rates.tvdt_ov.gt(0)].tvdt_ov.values * 12 / 365,
     growth_rates[growth_rates.tvdt_om.gt(0)].tvdt_om.values * 12 / 365]
plt.boxplot(X, labels=['Ovarian', 'Omental'])

# Add jitter
for e, i in enumerate(X):
    x = np.random.normal(e + 1, 0.04, size=len(i))
    plt.plot(x, i, 'k.', alpha=0.3)

plt.ylabel('TVDT (months)', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=14)
plt.yscale('log')
plt.tight_layout()
plt.savefig(path_for_output + 'figure_2_a.png', bbox_inches='tight')
plt.close()

'''
2. Box plots for omental and ovarian lesions by size
'''

# Get small and large lesions for ovarian and omental disease sites
om_vols = volumes[volumes.valid_om.ge(1)].groupby('anon_id').head(1)
ov_vols = volumes[volumes.valid_ov.ge(1)].groupby('anon_id').head(1)
tvdt_ov_small = growth_rates[growth_rates.anon_id.isin(ov_vols[ov_vols.vol_ov.le(ov_vols.vol_ov.median())].anon_id.unique())].tvdt_ov * 12 / 365
tvdt_ov_big = growth_rates[growth_rates.anon_id.isin(ov_vols[ov_vols.vol_ov.gt(ov_vols.vol_ov.median())].anon_id.unique())].tvdt_ov * 12 / 365
tvdt_om_small = growth_rates[growth_rates.anon_id.isin(om_vols[om_vols.vol_om.le(om_vols.vol_om.median())].anon_id.unique())].tvdt_om * 12 / 365
tvdt_om_big = growth_rates[growth_rates.anon_id.isin(om_vols[om_vols.vol_om.gt(om_vols.vol_om.median())].anon_id.unique())].tvdt_om * 12 / 365

# Creating boxplots
fig, ax = plt.subplots()
positions = [1, 2, 4, 5]
colors = ['blue', 'red']

bp1 = ax.boxplot([tvdt_ov_small, tvdt_om_small], positions=[1, 4], boxprops=dict(color=colors[0]))
bp2 = ax.boxplot([tvdt_ov_big, tvdt_om_big], positions=[2, 5], boxprops=dict(color=colors[1]))

# Add jitter
for e, i in enumerate([tvdt_ov_small, tvdt_ov_big, tvdt_om_small, tvdt_om_big]):
    x = np.random.normal(positions[e], 0.04, size=len(i))
    ax.plot(x, i, 'k.', alpha=0.3)


# Set x-ticks and x-tick labels
ax.set_xticks([1.5, 4.5])
ax.set_xticklabels(['Ovarian', 'Omental'], fontsize=14)

# Adding a legend
# blue_patch = plt.Line2D([0], [0], color='blue', lw=4)
# red_patch = plt.Line2D([0], [0], color='red', lw=4)
ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ['<= median', '> median'], loc='upper right', fontsize=16)
ax.set_ylabel('TVDT (months)', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=14)
ax.set_yscale('log')
plt.tight_layout()
plt.savefig(path_for_output + 'figure_2_b.png', bbox_inches='tight')
plt.close()

print('\n **** Stats TVDT small and big **** \n')
print('\n ******** OVARIAN - small ********** \n')
print(' {} \n'.format(tvdt_ov_small.agg('describe')))
print('\n ******** OVARIAN - big ********** \n')
print(' {} \n'.format(tvdt_ov_big.agg('describe')))
print('\n ******** OMENTAL - small ********** \n')
print(' {} \n'.format(tvdt_om_small.agg('describe')))
print('\n ******** OMENTAL - big ********** \n')
print(' {} \n'.format(tvdt_om_big.agg('describe')))

'''
p-value for the discrepancy between ovarian and omental TVDTs
'''
import scipy.stats as stats
tvdt_ov = growth_rates[growth_rates.tvdt_ov.ge(0)].tvdt_ov.values * 12 / 365
tvdt_om = growth_rates[growth_rates.tvdt_om.ge(0)].tvdt_om.values * 12 / 365


# Perform the Mann-Whitney U test - non-parametric, not normally distributed too little data
# First for om vs ov
statistic, p_value = stats.mannwhitneyu(tvdt_om, tvdt_ov, alternative='less')
print('**** OV vs OM *****')
print(f"U statistic: {statistic}")
print(f"p-value: {p_value}")

statistic, p_value = stats.mannwhitneyu(tvdt_ov_small.values, tvdt_ov_big.values, alternative='less')
print('**** OV big vs small *****')
print(f"U statistic: {statistic}")
print(f"p-value: {p_value}")

statistic, p_value = stats.mannwhitneyu(tvdt_om_small.values, tvdt_om_big.values, alternative='less')
print('**** OM big vs small *****')
print(f"U statistic: {statistic}")
print(f"p-value: {p_value}")

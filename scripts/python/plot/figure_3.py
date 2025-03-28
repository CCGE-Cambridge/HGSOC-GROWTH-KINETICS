'''
This is a script to plot the individual Gompertz curve fits for 11 patients with growing ovarian AND omental lesions along with their measured volumes.
We use a couple of them for figure 3 in the main manuscript

We also calculate the TTM and WOO mentioned in the Results section 'Median time to metastasis is 13 months for cases with growing primary and metastatic lesions'
'''
import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt


# Parameters
VMAX_OV = 5000 # cm3 from OV04 segmentations so far
VMAX_OM = 3000 # cm3 from OV04 segmentations so far
DETECTION_LIMIT_US = 0.5 # cm3 for reference
DETECTION_LIMIT_CA = 0.015 # cm3
V0 = 1e-9 # cm3
K_OV, K_OM = np.log(VMAX_OV/V0), np.log(VMAX_OM/V0) # Carrying capacity

# paths
path_to_gompertz_estimates = './../../../output/gompertz_params_{}.csv'
path_to_volumes = '../../../data/raw_volumes.csv'

path_for_output = '../../../output/plots/figure-3/'
if not os.path.exists(path_for_output):
    os.makedirs(path_for_output)

ov_params = pd.read_csv(path_to_gompertz_estimates.format('ov'), header=0, index_col=0)
om_params = pd.read_csv(path_to_gompertz_estimates.format('om'), header=0, index_col=0)
df_gompertz = pd.merge(ov_params, om_params, on='id', suffixes=('_ov','_om'), how='outer')
df_gompertz = df_gompertz.reset_index().rename(columns={df_gompertz.index.name:'anon_id'})
vols = pd.read_csv(path_to_volumes)

vols.sort_values(by=['anon_id', 'dt'], inplace=True)

'''
1. Plot individual Gompertz estimates for only those patients with valid ov and valid om (both increasing by 10% or more)
'''

df = df_gompertz.dropna()

def V(t, K, beta):
    return V0 * np.exp(K * (1 - np.exp(-beta * t)))

def get_time_to_vol_gompertz(V, beta, K, V0=1e-9):
    '''
    Assuming a Gompertz model, we estimate the time to reach a given volume, V

    Parameters
    ----------
    V: The volume for which we want the time. (cm3)
    beta: decay rate
    V0: init vol (cm3)
    K: log(Vinfty/V0)

    Returns
    t : time to reach the given volume --> t = -1/beta * log[ 1 - 1/K log[V/V0]]
    -------
    '''

    t = -1 / beta * np.log(1 - 1 / K * np.log(V / V0))

    return t

T = np.linspace(0, 10000, 1000)
for ix_, this_row in df.iterrows():
    V_ov = V(T, K_OV, this_row.loc['beta_ov'])
    V_om = V(T, K_OM, this_row.loc['beta_om'])

    pat_vols = vols[vols.anon_id == this_row.loc['anon_id']]
    t_ = [this_row.loc['t1_ov'] , this_row.loc['t1_ov'] + pat_vols.dt.iloc[-1]]
    t_ = [i * 12 / 365 for i in t_]
    dt = this_row.loc['t1_ov'] - this_row.loc['t1_om'] # Gap in time beteween initiation of ovarian and omental tumours

    # Creating the figure and primary axis
    fig, ax1 = plt.subplots()

    # Plotting ovarian tumour volumes on the primary y-axis
    color = 'tab:blue'
    ax1.set_xlabel('Time (months)', fontsize=18)
    ax1.set_ylabel('Ovarian burden ($cm^3$)', color=color, fontsize=18)
    ax1.plot(T * 12 / 365, V_ov, color=color, label='ov', linewidth=2)
    ax1.plot(t_, 1e-3 * pat_vols['vol_ov'], linestyle='None', marker = 'x', ms = 12, mec = color)
    
    # Plot reference detection limit
    ax1.plot([-100, T[-1] * 12 / 365], [DETECTION_LIMIT_US, DETECTION_LIMIT_US], 'k-')
    ax1.plot([-100, T[-1] * 12 / 365], [DETECTION_LIMIT_CA, DETECTION_LIMIT_CA], color='gray')

    ax1.tick_params(axis='y', labelcolor=color, labelsize=14)
    ax1.tick_params(axis='x', labelsize=14)
    ax1.set_ylim([1e-9, (pat_vols['vol_om'] + pat_vols['vol_ov']).max() * 1.5 * 1e-3])
    ax1.set_yscale('log')

    # Get handles and labels from ax1
    handles1, labels1 = ax1.get_legend_handles_labels()
    # Create the legend for ax1 outside the plot, on the top
    ax1.legend(handles1, labels1, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=len(labels1), prop={'size': 14})

    # Plotting omental volumes on the secondary y axis
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Omental burden ($cm^3$)', color=color, fontsize=18)
    ax2.plot((T + dt) * 12 / 365, V_om, color=color, label='om', linewidth=2)
    ax2.plot(t_, 1e-3 * pat_vols['vol_om'], linestyle='None', marker = 'x', ms = 12, mec = color)
    ax2.tick_params(axis='y', labelcolor=color, labelsize=14)
    ax2.tick_params(axis='x', labelsize=14)
    ax2.set_ylim([1e-9, (pat_vols['vol_om'] + pat_vols['vol_ov']).max() * 1.5 * 1e-3])
    ax2.set_yscale('log')

    fig.tight_layout()
    plt.xlim(right=(t_[-1] * 1.1))
    plt.xlim(left=min([0 , (T[0] + dt) * 12 / 365]))
    plt.savefig(path_for_output + '{}.png'.format(this_row.loc['anon_id']))
    plt.close()

'''
Now calculate 
1. Time between two tumours (TTM)
2. WOO for detection
'''

# 1. diff between t1_ov and t1_om for all cases with valid ovarian AND omental lesions
df_both = df.copy()
df_both.loc[:,'d_t1'] = df_both.t1_ov - df_both.t1_om
df_both.loc[:,'d_t1_abs'] = np.abs(df_both.t1_ov - df_both.t1_om)

print('Stats for time between primary and secondary sites (months) \n')
print(df_both.d_t1_abs.apply(lambda x: x * 12 / 365).agg('describe'))

# 2. WOO between the primary lesion reaching the detection limit and the onset of metastasis
df_both.loc[:,'t_detect_ov_us'] = df_both.apply(lambda x: get_time_to_vol_gompertz(DETECTION_LIMIT_US,
                                                                                 x.loc['beta_ov'],
                                                                                 K_OV,
                                                                                 V0),
                                             axis=1)
df_both.loc[:,'t_detect_om_us'] = df_both.apply(lambda x: get_time_to_vol_gompertz(DETECTION_LIMIT_US,
                                                                                x.loc['beta_om'],
                                                                                K_OM,
                                                                                V0),
                                             axis=1)
df_both.loc[:,'t_detect_ov_ca'] = df_both.apply(lambda x: get_time_to_vol_gompertz(DETECTION_LIMIT_CA,
                                                                                 x.loc['beta_ov'],
                                                                                 K_OV,
                                                                                 V0),
                                             axis=1)
df_both.loc[:,'t_detect_om_ca'] = df_both.apply(lambda x: get_time_to_vol_gompertz(DETECTION_LIMIT_CA,
                                                                                x.loc['beta_om'],
                                                                                K_OM,
                                                                                V0),
                                             axis=1)

# Note: the WOO isn't t_detect_ov - t1_om because t1_ov/om is time since ovarian/omental lesion started to the first scan
# It is instead t1_ov - t_detect_ov_us - t1_om
# This is easier to understand when looking at figure 3.
df_both.loc[:, 'WOO_met_init_us'] = df_both.apply(lambda x: np.max([x['t1_ov'] - x['t_detect_ov_us'] - x['t1_om'],
                                                            x['t1_om'] - x['t_detect_om_us'] - x['t1_ov']]),
                                          axis=1)
df_both.loc[:, 'WOO_met_init_ca'] = df_both.apply(lambda x: np.max([x['t1_ov'] - x['t_detect_ov_ca'] - x['t1_om'],
                                                            x['t1_om'] - x['t_detect_om_ca'] - x['t1_ov']]),
                                          axis=1)
df_can_detect_us = df_both[df_both.WOO_met_init_us.gt(0)]
df_can_detect_ca = df_both[df_both.WOO_met_init_ca.gt(0)]

print('\n Stats for WOO between primary detection (ULTRASOUND) and metastasis (months) \n')
print(df_can_detect_us.WOO_met_init_us.apply(lambda x: x * 12 / 365).agg('describe'))
print('\n Stats for WOO between primary detection (CA125) and metastasis (months) \n')
print(df_can_detect_ca.WOO_met_init_ca.apply(lambda x: x * 12 / 365).agg('describe'))

'''
Get the size of the primary at the time the secondary is initiated
This vector will be used in simulations to draw the time at which metastasis will occur for the different lesions.
'''
df_size = df.copy()
df_size.loc[:,'d_t1'] = df_size.t1_ov - df_size.t1_om
df_size['size_at_met'] = 0
def size_at_met(row):
    if row.d_t1 > 0:
        return V(row.loc['d_t1'], K_OV, row.loc['beta_ov'])
    else:
        return V(-1 * row.loc['d_t1'], K_OM, row.loc['beta_om'])

df_size.loc[:,'size_at_met'] = df_size.apply(lambda x: size_at_met(x), axis=1)

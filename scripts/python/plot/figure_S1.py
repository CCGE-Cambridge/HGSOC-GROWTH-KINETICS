'''
Script to just plot all the volumes from the 37 patients..
the ones we discard will be in red lines and red crosses/dots
'''
import pandas as pd
import matplotlib.pyplot as plt

# Paths
path_to_volumes = './../../../data/raw_volumes.csv'
path_for_output = './../../../output/plots/suppmat/'

vols = pd.read_csv(path_to_volumes)
vols = vols.sort_values(by=['anon_id', 'dt'])

# Plotting for ovarian lesions
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
for pat_id, pat_data in vols.groupby('anon_id'):

    # Get time
    dt = pat_data.dt.iloc[-1]
    T = [0, (dt * 12 / 365)]

    # Plot ovarian
    if pat_data.iloc[-1].valid_ov:
        if pat_data.iloc[-1].valid_om:
            ax1.plot(T, pat_data['vol_ov'] * 1e-3, marker='o', color='blue')
        else:
            ax1.plot(T, pat_data['vol_ov'] * 1e-3, marker='o', color='black')
    else:
        ax1.plot(T, pat_data['vol_ov'] * 1e-3, marker='o', color='red')

    # Plot omental
    if pat_data.iloc[-1].valid_om:
        if pat_data.iloc[-1].valid_ov:
            ax2.plot(T, pat_data['vol_om'] * 1e-3, marker='o', color='blue')
        else:
            ax2.plot(T, pat_data['vol_om'] * 1e-3, marker='o', color='black')
    else:
        ax2.plot(T, pat_data['vol_om'] * 1e-3, marker='o', color='red')

for ax in [ax1, ax2]:
    ax.set_xlabel('Time (months)', fontsize=18)
    ax.set_ylabel('Volume ($cm^3$)', fontsize = 18)
    ax.tick_params(axis='x', labelsize=16)  # Change xtick label size
    ax.tick_params(axis='y', labelsize=16)  # Change xtick label size
    ax.set_yscale('log')
    # ax.set_ylim([ax.get_ylim()[0], 3e3])

fig1.tight_layout()
fig2.tight_layout()
fig1.savefig(path_for_output + 'figure_S1_a.png')  # Save the figure
fig2.savefig(path_for_output + 'figure_S1_b.png')  # Save the figure
plt.close()

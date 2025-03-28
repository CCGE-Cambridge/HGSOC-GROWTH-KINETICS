"""
Script to plot violin plots for all the vmax sensitivity runs from Matlab
"""
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


Vmaxs = np.linspace(2000, 10000, 17)

# Paths
path_to_sensitivity = './../../../output/sensitivity-analysis/vmax_{}.csv'
path_for_output = './../../../output/plots/suppmat/'
if not os.path.exists(path_for_output):
    os.makedirs(path_for_output)

sites = {'ov': 5000,
         'om': 3000}
fig_names = {'ov': 'c',
         'om': 'd'}
for site, vmax in sites.items():
    df = pd.read_csv(path_to_sensitivity.format(site))
    data = [np.random.lognormal(mean, sigma, 50000) for mean, sigma in zip(df['mean-logt1'].values, df['std-logt1'].values)]

    # Creating a DataFrame for plotting
    df_plot = pd.DataFrame({f'{int(Vmaxs[i])}': dist for i, dist in enumerate(data)})

    plt.violinplot(df_plot * 12 / 365, showmedians=True, showextrema=False, points=500) #, quantiles=[[0.25, 0.75]]*20)
    plt.ylim([0, 100])
    plt.xticks(ticks=np.arange(1, len(Vmaxs)+1), labels=['{}'.format(int(i)) for i in Vmaxs], rotation=45, fontsize=16)

    # Plot a horizontal line at the median value of either 5000 or 3000
    y_val = df_plot[str(vmax)].median() * 12 / 365
    plt.plot([0.5, 17.5], [y_val, y_val], 'r--', linewidth=1)

    # Customizing the plot
    plt.xlabel('$V_{\infty}$', fontsize=18)
    plt.ylabel('$t_1$ (months)', fontsize=18)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(path_for_output + '/figure_S2_{}.png'.format(fig_names[site]), bbox_inches='tight')
    plt.close()

    # Get stats on the coefficient of variation of mean-logt1 and std-logt1
    print('CV for {} for the mean of log t1 is {}'.format(site, np.std(np.exp(df['mean-logt1'])) / np.mean(np.exp(df['mean-logt1']))))

    # Get stats on the coefficient of variation of mean-logt1 and std-logt1 for the choices above 5000 and 3000 cm3 for ovarian and omental
    indices = np.where(Vmaxs >= vmax)[0]
    df_ = df.iloc[indices]
    print('CV for {} with vmax >= {} for the mean of log t1 is {}'.format(site, vmax, np.std(np.exp(df_['mean-logt1'])) / np.mean(np.exp(df_['mean-logt1']))))
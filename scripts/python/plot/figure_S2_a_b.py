"""
Script to plot violin plots for all the measurement error sensitivity runs from Matlab
"""
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# Paths
path_to_sensitivity = './../../../output/sensitivity-analysis/measurement_10per_{}.csv'
path_for_output = './../../../output/plots/suppmat'
if not os.path.exists(path_for_output):
    os.makedirs(path_for_output)

# data for lines corresponding to estimated fixed effects of t1 using the actual data
sites = {'ov': 20.4024,
         'om': 14.6871}
figure_labels = {'ov': 'a',
                 'om': 'b'}

# Go through each site and each run and plot the distribution of t1 after using the log(t1) values
for site, val in sites.items():
    df = pd.read_csv(path_to_sensitivity.format(site))

    # Check if parameters of log(t1) estimated by actually drawing from a distribution of log(q) and log(beta)
    # and getting distribution of log(t1) as log(q) - log(beta) across 100,000 iterations in Matlab
    # are the same as when derived from params of log(beta) and log(q) from Matlab
    # https://online.stat.psu.edu/stat500/book/export/html/572
    # Yes they are! All are less than 1% error difference
    df['mu_est'] = df['fe-logq'] - df['fe-logbeta']
    df['diff_mu'] = np.abs(df['mu_est'] - df['mean-logt1']) / df['mean-logt1']
    df['std_est'] = np.sqrt(df['sd-logq'] + df['sd-logbeta'])
    df['diff_std'] = np.abs(df['std_est'] - df['std-logt1']) / df['std-logt1']

    # Now draw from the distribution of log-t1 with 50,000 samples
    data = [np.random.lognormal(mean, sigma, 50000) for mean, sigma in zip(df['mean-logt1'].values, df['std-logt1'].values)]

    # Creating a DataFrame for plotting
    df_plot = pd.DataFrame({f'Dist {i + 1}': dist for i, dist in enumerate(data)})

    # Create a violin plot
    plt.violinplot(df_plot * 12 / 365, showmedians=True, showextrema=False, points=500) #, quantiles=[[0.25, 0.75]]*20)
    plt.ylim([0, 100])
    plt.xticks([])

    # Plot the median from the no-noise case taken directly from the data
    plt.plot([0.5, 20.5], [val, val], 'r--', linewidth=1)

    # Customizing the plot
    plt.xlabel('Runs', fontsize=18)
    plt.ylabel('$t_1$ (months)', fontsize=18)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(path_for_output + '/figure_S2_{}.png'.format(figure_labels[site]), bbox_inches='tight')
    plt.close()

    # Get stats on the coefficient of variation of mean-logt1 and std-logt1
    print('CV for {} for the mean of log t1 is {}'.format(site, np.std(np.exp(df['mean-logt1'])) / np.mean(np.exp(df['mean-logt1']))))
    # print('CV for {} for the std of log t1 is {}'.format(site, np.std(df['std-logt1']) / np.mean(df['std-logt1'])))


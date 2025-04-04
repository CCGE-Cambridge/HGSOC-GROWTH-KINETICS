'''
This is a script to simulate a population using the Gompertz growth parameters
calculated using the NLME model.
I will also use the size of the primary tumour at which metastasis was estimated to have
started from the data (11 patients with growing lesions in both sites)
'''
import pandas as pd
import numpy as np
import os
import sys
import random


# Paths
path_for_output = './../../output/simulations/'
if not os.path.exists(path_for_output):
    os.makedirs(path_for_output)

'''
Parameter setting
'''
N_TUMOURS = 10000 # This is thue number of tumours
V0 = 1e-9 # starting volume
LIMIT_US = 0.5 # ultrasound detection limit (cm3)
LIMIT_CA = 0.015 # ca125 detection limit (cm3)

# Max volumes set before hand while doing NLME
# mean and stds of log(beta) estimated using NLME - from Matlab
# beta is the Gompertz decay rate
ovarian_params = {'vmax': 5000, # max volume -> determines K = log(Vmax/V0) where V0 = 1e-9
                  'ln_beta_mean': -5.7905, # mean of log(beta) from the NLME model
                  'ln_beta_std': np.sqrt(0.8802),} # std of log(beta) from the NLME model. It actually gives us the variance.
omental_params = {'vmax': 3000,
                  'ln_beta_mean': -5.5188,
                  'ln_beta_std': np.sqrt(0.8799),}

# sizes at metastasis from the 11 cases with growing lesions in both sites (got this from plot_gompertz_indi_from_matlab.py)
pt_size_at_met = [6.58634689e-07, 7.49423904e-04, 1.39645166e-02, 8.36070562e-02,
       3.88565349e-01, 1.04536906e+00, 1.22672942e+00, 2.45892441e+00,
       3.67071643e+00, 3.87580839e+00, 3.99019334e+00]
mets = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] # indices

# Probability of omental primary vs ovarian primary - 4 out of 11 had omental disease first
p_om = 4/11

# Gompertz function
def V(t, K, beta):
    '''

    Parameters
    ----------
    t : time
    K : carrying capacity parameter, K = ln(Vmax/V0)
    beta : decay rate

    Returns
    -------
    V = V0 * exp(K * (1 - exp(-beta * t))
    '''
    return V0 * np.exp(K * (1 - np.exp(-beta * t)))


def get_time_to_vol_gompertz(V, beta, K, V0=1e-9):
    '''
    Assuming a Gompertz model, we estimate the time to reach a given volume, V

    Parameters
    ----------
    V: The volume for which we want the time. (cm3)
    beta: decay rate
    V0: init vol (cm3)
    K: log(Vmax/V0)

    Returns
    t : time to reach the given volume --> t = -1/beta * log[ 1 - 1/K log[V/V0]]
    -------
    '''

    t = -1 / beta * np.log(1 - 1 / K * np.log(V / V0))

    return t

list_dicts = []
# Loop through each tumour simulation
for i in range(N_TUMOURS):

    # Uniform probability used to determine if it's an ovarian/omental primary
    omental = False
    if np.random.uniform() <= p_om:
      omental=True
      params_pt = omental_params
      params_met = ovarian_params
    else:
      params_pt = ovarian_params
      params_met = omental_params

    # Draw the size of the primary tumour (PT) at the onset of metastasis from the estimated size of PT at metastasis
    # for the 11 cases with growing omental and ovarian lesions
    size_at_met = random.choice(pt_size_at_met)

    # Draw the primary and metastatic decay rates from their respective distributions
    # When the ovarian is the primary, omental is the metastatic site and vice versa
    beta_pt = np.exp(np.random.normal(params_pt['ln_beta_mean'], params_pt['ln_beta_std']))
    beta_met = np.exp(np.random.normal(params_met['ln_beta_mean'], params_met['ln_beta_std']))

    # Set the K parameter
    K = np.log(params_pt['vmax'] / V0)

    # Time taken for PT to metastasise
    t_to_met = get_time_to_vol_gompertz(size_at_met, beta_pt, K, V0)

    # Time taken for PT to reach the US / CA125 detection limit
    t_to_detect_CA125 = get_time_to_vol_gompertz(LIMIT_CA, beta_pt, K, V0)
    t_to_detect_US = get_time_to_vol_gompertz(LIMIT_US, beta_pt, K, V0)

    # Size of mets at CA125 / US detection limit
    K_met = np.log(params_met['vmax'] / V0)
    met_size_at_CA_detect = V(t_to_detect_CA125 - t_to_met, K_met, beta_met)
    met_size_at_US_detect = V(t_to_detect_US - t_to_met, K_met, beta_met)

    # store information
    dict_ = {'ix': i, # index of current tumour simulation
             'omental': omental, # whether it is the omental primary
             'size_at_met': size_at_met,
             'beta_pt': beta_pt,
             'beta_met': beta_met,
             'time_to_met': t_to_met,
             'time_to_ca125': t_to_detect_CA125,
             'time_to_US': t_to_detect_US,
             'met_size_at_ca125': met_size_at_CA_detect,
             'met_size_at_US': met_size_at_US_detect,
             }
    list_dicts.append(dict_)

df = pd.DataFrame(list_dicts)
df.to_csv(path_for_output + 'sims.csv')

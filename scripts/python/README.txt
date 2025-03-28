This folder contains the python scripts used in this paper.

1. I_calculate_tvdts.py
-- This script takes in the raw volumes in the data folder and calculated the tumour volume doubling times

2. II_simulate_population.py
-- This script simulates 10,000 tumours by drawing from the distribution of the Gompertz parameters for ovarian and omental lesions.
-- It relies on the results obtained in the MATLAB section but can be run directly.

3. III_analyse_simulation_results.py
-- This script takes in the results from II_simulate_population.py
-- It outputs stats on the number of tumours that will metastasise before detection and the WOO for US/CA125 based detection for
   those tumours that can be detected before metastasis.
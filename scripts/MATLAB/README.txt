This is the folder with all the MATLAB scripts used for the analyses.

1. run_fitting_gompertz.m:
-- This is the script for fitting a nonlinear mixed effects model with a Gompertz growth assumption to the data
-- The script takes in the raw volumes in the data folder
-- It outputs the individual parameter estimates for the 24 cases with growing ovarian lesions
-- It also provides the fixed effects parameters (phi) and the standard deviation of the random effects (PSI) for
   the decay rate beta (phi(1) and PSI(1,1)) and the product q= beta*t1 (phi(2) and PSI(2,2))
-- We also calculate the confidence interval for the exponential of the geometric mean (log(t1^pop))
-- The results of this are used for figure 3 and for the results in the section 'Median time to metastasis is 13 months for cases with growing primary and metastatic lesions'

2. run_fitting_exponential.m:
-- This is the script for fitting a nonlinear mixed effects model with an exponential growth assumption to the data
-- The script takes in the raw volumes in the data folder
-- It outputs the individual parameter estimates for the 24 cases with growing ovarian lesions
-- It also provides the fixed effects parameters (phi) and the standard deviation of the random effects (PSI) for 
   the growth rate mu (phi(1) and PSI(1,1)) and the product of mu and t1 (phi(2) and PSI(2,2))
-- We also calculate the confidence interval for the exponential of the geometric mean (log(t1^pop))

3. run_measurement_sensitivity.m:
-- This script applies 10% randomly distributed Gaussian noise to the data and does the fitting
-- This is done 20 times
-- It returns a 20x6 vector with the means and standard deviations of log(t1), log(beta), and log(q=beta*t1).
   The latter two are estimated directly by the NLME while the log(t1) is estimated by generating 100,000 values of log(beta)
   and log(q) from their fixed effects values and the SD of their random effects, and then calculating the mean and SD of log(t1) = log(q) - log(beta)
-- The results of this are used to plot figure S2 a and b

4. run_vmax_sensitivity.m
-- This script is used to study the changes in t_1 as the Vmax value is changed from 2,000 to 10,000 cm3
-- It returns a 20x6 vector with the means and standard deviations of log(t1), log(beta), and log(q=beta*t1).
   The latter two are estimated directly by the NLME while the log(t1) is estimated by generating 100,000 values of log(beta)
   and log(q) from their fixed effects values and the SD of their random effects, and then calculating the mean and SD of log(t1) = log(q) - log(beta)
-- The results of this are used to plot figure S2 c and d
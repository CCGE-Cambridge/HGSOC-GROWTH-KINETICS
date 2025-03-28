clear all
close all
clc 

% % This is a script to fit a gompertz model using NLME to OV04 data
% Parameters and model
% use 'om' for omental and 'ov' for pelvic/ovarian disease
data.TYPE = 'om';
data.V0 = 1e-9;

% Load data
T_raw = readtable(['raw_volumes.csv']);
if data.TYPE == 'ov'
    T = T_raw(T_raw.valid_ov == 1, :);
    [data.ids, data.volumes, data.dt] = deal(T.anon_id, T.vol_ov * 1e-3, T.dt);
    data.Vmax = 5000;
else
    T = T_raw(T_raw.valid_om == 1, :);
    [data.ids, data.volumes, data.dt] = deal(T.anon_id, T.vol_om * 1e-3, T.dt);
    data.Vmax = 3000;
end

% Fit the model
[phi, PSI, stats, br, y] = fit_gompertz(data);

% Get and store statistics for the model
% [mean_t1, std_t1] = get_stats(phi, PSI, 0, data);

% save(['./output/general/results_' data.TYPE])

% Output just the individual parameters for use in python
out.id = data.ids(1:2:end);
out.beta = exp(phi(1)+ br(1,:))';
out.t1 = exp(phi(2) + br(2,:))' ./ out.beta(:);
writetable(struct2table(out), ['./output/gompertz_params_' data.TYPE '.csv'])

%% CI for exp(ln(t1_pop)) from fixed effects estimates of log(beta) and log(beta * t1)
% Derivation in the supplementary material
% t statistic threshold from tinv based on N_obs
n_obs = size(br, 2);
mu_ln_t1 = phi(2) - phi(1); % pop level log(t1) = log(beta_fe*t1_fe) - log(beta_fe)
SE_ln_t1 = sqrt(stats.sebeta(1)^2 + stats.sebeta(2)^2); % standard error of difference of normally distributed variables
t_stat = tinv(0.975, n_obs-1); % t statistic corresponding to the 95% CI
CI_ln_t1 = [mu_ln_t1 - t_stat*SE_ln_t1, mu_ln_t1 + t_stat*SE_ln_t1];
CI_t1 = exp(CI_ln_t1) * 12 / 365; % convert to months
FE_t1 = exp(phi(2)) / exp(phi(1)) * 12 / 365;
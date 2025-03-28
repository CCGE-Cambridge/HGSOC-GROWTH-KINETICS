%% Script to run volume sensitivity to Vmax
% This is a script that varies the Vmax and fits the Gompertz NLME model for each choice of Vmax
clc
clear all

% Parameters and model
data.TYPE = 'ov';
data.V0 = 1e-9;
Vmaxs = 2000:500:10000;

% Load data
T_raw = readtable(['raw_volumes.csv']);
if data.TYPE == 'ov'
    T = T_raw(find(T_raw.valid_ov == 1), :);
    [data.ids, data.volumes, data.dt] = deal(T.anon_id, T.vol_ov * 1e-3, T.dt);
    data.Vmax = 5000;
else
    T = T_raw(find(T_raw.valid_om == 1), :);
    [data.ids, data.volumes, data.dt] = deal(T.anon_id, T.vol_om * 1e-3, T.dt);
    data.Vmax = 3000;
end

means = [];
stds = [];
febeta = [];
feq = [];
sdbeta = [];
sdq = [];

% Loop through different Vmaxs
for i = 1:length(Vmaxs)
    data.Vmax = Vmaxs(i);
    
    % Fit the model
    [phi, PSI, stats, br] = fit_gompertz(data);
    close();

    % Get and store mean and std of log(t1) from the generated distribution
    [mean_t1, std_t1] = get_stats(phi, PSI, 0, data);
    means(i) = mean_t1;
    stds(i) = std_t1;

    % Store the mean and std of log(q) and log(beta) - will be used to
    % estimate the mean and std of log(t1) from the formula for a
    % difference of two normally distributed variables
    febeta(i) = phi(1);
    feq(i) = phi(2);
    sdbeta(i) = PSI(1,1);
    sdq(i) = PSI(2,2);
end

% Create a table with the results
resultsTable = table(means', stds', febeta', feq', sdbeta', sdq' ,'VariableNames', {'mean-logt1', 'std-logt1', 'fe-logbeta', 'fe-logq', 'sd-logbeta', 'sd-logq'});

% Write the table to a CSV file
writetable(resultsTable, ['./output/sensitivity-analysis/vmax_' data.TYPE '.csv']);


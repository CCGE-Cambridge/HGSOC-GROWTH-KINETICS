%% Sensitivity to measurement errors
% This is a script that adds 10% random Gaussian noise to the measured
% volumes and calculates the fixed effects ln(t1_pop). This is run 20 times
% to see how the estimates population geometric mean varies.
clear all
clc 

% Parameters
data.TYPE = 'om';
data.V0 = 1e-9;
percent_noise = 10;
n_evals = 20;

% Load data
T_raw = readtable(['raw_volumes.csv']);
if data.TYPE == 'ov'
    T = T_raw(find(T_raw.valid_ov == 1), :);
    [data.ids, volumes, data.dt] = deal(T.anon_id, T.vol_ov * 1e-3, T.dt);
    data.Vmax = 5000;
else
    T = T_raw(find(T_raw.valid_om == 1), :);
    [data.ids, volumes, data.dt] = deal(T.anon_id, T.vol_om * 1e-3, T.dt);
    data.Vmax = 3000;
end

means = [];
stds = [];
febeta = [];
feq = [];
sdbeta = [];
sdq = [];
for i = 1:n_evals
    disp(i);
    flag = 0;
    while flag < 1        
        noise = normrnd(0, percent_noise/100, [length(volumes),1]);
        data.volumes = volumes .* (1 + noise);
        
        % Check if all have an increase! Otherwise, this is a biased result as
        % they would not be included in the analysis (we would only include
        % those with an increase in volumes)
        diff = data.volumes(2:2:end) - data.volumes(1:2:end);
        if min(diff) > 0
            flag = 1;
        end 
    end 

    % Fit the model
    [phi, PSI, stats, br] = fit_gompertz(data);
    close();

    % Get and store mean and std of log(t1)
    [mean_t1, std_t1] = get_stats(phi, PSI, 0, data);
    means(i) = mean_t1;
    stds(i) = std_t1;
    febeta(i) = phi(1);
    feq(i) = phi(2);
    sdbeta(i) = PSI(1,1);
    sdq(i) = PSI(2,2);
end

% Create a table with the results
resultsTable = table(means', stds', febeta', feq', sdbeta', sdq' ,'VariableNames', {'mean-logt1', 'std-logt1', 'fe-logbeta', 'fe-logq', 'sd-logbeta', 'sd-logq'});

% Write the table to a CSV file
writetable(resultsTable, ['./output/sensitivity-analysis/measurement_' num2str(percent_noise) 'per_' data.TYPE '.csv']);

% Get the stats of variation / median
CV_t1 = std(exp(means)) / mean(exp(means));
CV_std_t1 = std(stds) / mean(stds);
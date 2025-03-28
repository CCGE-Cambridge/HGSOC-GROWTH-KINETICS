function [mean_t1, std_t1] = get_stats(phi, PSI, PLOT_HIST, data)
% GET_STATS function to get the mean and std of the log normal distributions
% We sample from distrbutions of beta and beta*t1 and then calculate t1
% we then get the mean and std of log(t1) because this is normally
% distributed

[fe_beta, fe_q] = deal(phi(1), phi(2));
[sd_beta, sd_q] = deal(sqrt(PSI(1,1)), sqrt(PSI(2,2))); % q here is -beta * t1

% Both beta and q are log transformed so 
beta_vec = exp(normrnd(fe_beta, sd_beta ,[1, 100000]));
q_vec = exp(normrnd(fe_q, sd_q ,[1, 100000])); 
t1_vec = q_vec ./ beta_vec; % Here q = beta * t1

log_t1_vec = log(t1_vec);
mean_t1 = mean(log_t1_vec);
std_t1 = std(log_t1_vec);

% Plot a histogram of t1
if PLOT_HIST == 1
    t1_vec = t1_vec * 12 / 365;
    histogram(t1_vec); hold on;
    yLimits = ylim;
    plot([median(t1_vec), median(t1_vec)], yLimits, 'r--', 'LineWidth', 2);
    xlim([0 10000 * 12 / 365]);
    title(['Median t1 =' num2str(median(t1_vec)) 'months'])
    xlabel('t1 (Months)', 'FontSize', 14)
    saveas(gcf, ['./output/figures/hist_t1_' data.TYPE '.fig'])
    saveas(gcf, ['./output/figures/hist_t1_' data.TYPE '.png'])
end

end


function [phi_final, PSI_final, stats_final, br_final, y] = fit_gompertz(data, full_cov)
    % This is a function to fit a gompertz model using NLME to OV04 data
    % Assumptions:
    % Use the formulation: ln (1 - omega/K) = ln(q) - beta * dt
    %                       q: beta * t1 
    %                       beta: decay rate
    %                       t1: time between lesion initiation and the first scan
    %
    % ------ INPUTS -----
    % data: object that stores
    %       Vmax: maximum disease volume (cm3)
    %       V0 : Volume at time 0 = 1e-9 (cm3)
    %       volumes : measured volumes (cm3)
    %       dt : time interval between volume measurements (days)
    % full_cov: whether we want full covariance or not (1 or 0)
    % log_param: whether we want to use log transformed parameters (1 or 0)
    %
    % ------ OUTPUTS ------
    % phi_final: 2x1 array that contains the estimated fixed effects of
    %           log(beta) and log(beta * t1)
    % PSI_final: 2x2 covariance matrix of the random effects
    % stats_final: summary of the statistics including the errors
    % br_final: 2 x N array where N is the number of patients with the
    %           individual difference in log(beta) and log(beta * t1) with
    %           respect to the fixed effects
    % y : the transformed observations, where y = ln (1 - omega/K) and
    %     omega = log(V / V0)
    
    if nargin<2
        full_cov = 0;
    end
    
    K = log(data.Vmax / data.V0);
    
    %% Transformed Gompertz model
    % phi(1) is beta
    % phi(2) = -ln (q) where q = beta * t1 where t1 is the time to
    % We have -ln(1-omega/K) = -ln(q) + beta * dt
    % Minus sign on both sides to force ln (q) to be +ve
    omega = log(data.volumes./data.V0);
    y = -1*log(1 - omega./K); 
    model = @(phi, dt)(phi(2) + phi(1).*dt);
    
    % Starting condition - note that since we are log transforming the
    % coordinates, we should also log transform the starting guess
    % beta is around 0.005 and t1 is around 200 days
    beta_0 = 0.005;
    t1_0 = 500;
    q_0 = beta_0 * t1_0; 
    
    phi0 = [log(beta_0); 
                log(q_0)];

    xform = [1, 1]; % log transform parameters
    
    if full_cov == 1
        cov_pattern = ones(2);
    else
        cov_pattern = eye(2);
    end

    %% NLME (linearized)
    options = statset('MaxIter', 1000, 'TolFun', 1e-6, 'TolX', 1e-6);
    
    [phi_final, PSI_final, stats_final, br_final] = nlmefit(data.dt ,y ,data.ids,[], model, phi0, ...
        'RefineBeta0','off','ParamTransform',xform, 'Options', options, ...
        'CovPattern',cov_pattern, ...
        'ErrorModel','constant'); 

end
   

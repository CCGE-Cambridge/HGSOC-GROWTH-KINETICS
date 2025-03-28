function [phi_final, PSI_final, stats_final, br_final] = fit_exponential(data)
    % This is a function to fit an exponential model using NLME to OV04 data
    % Assumptions:
    % 1. V0 = 1e-9 like for Gompertz
    % Use the formulation: ln (V/V0) = mu * (t1 + dt)
    % ------ INPUTS ------
    % data: object that stores
    %       Vmax: maximum disease volume (cm3)
    %       V0 : Volume at time 0 = 1e-9 (cm3)
    %       volumes : measured volumes (cm3)
    %       dt : time interval between volume measurements (days)
    %
    % ------ OUTPUTS ------
    % phi_final: 2x1 array that contains the estimated fixed effects of
    %           log(beta) and log(beta * t1)
    % PSI_final: 2x2 covariance matrix of the random effects
    % stats_final: summary of the statistics including the errors
    % br_final: 2 x N array where N is the number of patients with the
    %           individual difference in log(beta) and log(beta * t1) with
    %           respect to the fixed effects

      
    %% Transformed problem
    % phi(1) is mu
    % phi(2) = mu * t1
    % We have ln (V/V0) = phi(2) + phi(1) * dt
    omega = log(data.volumes./data.V0);
    y = omega; 
    model = @(phi, dt)(phi(2) + phi(1).*dt);
    

    % Starting condition - note that since we are log transforming the
    % coordinates, we should also log transform the starting guess
    % growth rate mu is around 0.02 and t1 is around 500 days
    mu_0 = 0.02;
    t1_0 = 500;
    q_0 = mu_0 * t1_0; % -ln(q) = beta_0 * t1
    
    phi0 = [log(mu_0); 
        log(q_0)];
    xform = [1, 1]; % log transform parameters to ensure only positive values
    
    %% NLME (linearized)
    [phi_final, PSI_final, stats_final, br_final] = nlmefit(data.dt ,y ,data.ids,[], model, phi0, ...
        'RefineBeta0','off','ParamTransform',xform, 'Options', statset('MaxIter', 5000));
end
   
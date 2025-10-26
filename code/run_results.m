%% run_results.m  — simulate, create obs, and plot

% --- Parameters (COVID-like) ---
N = 1e6; beta = 0.39; sigma = 0.25; gamma = 0.17; alpha = 0.002;
Tend = 180;
y0_seir  = [N-10; 5; 5; 0];
y0_seird = [N-10; 5; 5; 0; 0];

% --- SEIR ODE ---
SEIR = @(t,y)[ ...
   -beta*y(1)*y(3)/N; ...
    beta*y(1)*y(3)/N - sigma*y(2); ...
    sigma*y(2) - gamma*y(3); ...
    gamma*y(3)];

% --- SEIRD ODE ---
SEIRD = @(t,y)[ ...
   -beta*y(1)*y(3)/N; ...
    beta*y(1)*y(3)/N - sigma*y(2); ...
    sigma*y(2) - (gamma+alpha)*y(3); ...
    gamma*y(3); ...
    alpha*y(3)];

% --- Integrate both (ode45) ---
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t_seir,  y_seir]  = ode45(@(t,y)SEIR(t,y),  [0 Tend], y0_seir,  opts);
[t_seird, y_seird] = ode45(@(t,y)SEIRD(t,y), [0 Tend], y0_seird, opts);

% --- Build synthetic observed incidence from SEIRD + Poisson noise ---
t_daily = (0:Tend)';                               % daily grid
E_truth  = interp1(t_seird, y_seird(:,2), t_daily, 'pchip');  % exposed
inc_true = sigma * E_truth;                        % model incidence (lambda)
rng(1);                                            % reproducible
obs = mypoissrnd(max(inc_true,0));                 % daily counts (column)
obs = obs(:);

% --- Call graphs() to plot + metrics ---
[rmseS, rmseD, aicS, aicD] = graphs(obs, t_seir, y_seir, t_seird, y_seird, sigma);

% --- Display metrics in Command Window ---
fprintf('RMSE  SEIR : %.2f\nRMSE  SEIRD: %.2f\nAIC   SEIR : %.1f\nAIC   SEIRD: %.1f\n', ...
        rmseS, rmseD, aicS, aicD);
% After graphs(...)
plot_cumulative(obs, t_seir, y_seir, t_seird, y_seird, sigma, Tend);
Talpha = sensitivity_alpha(obs, beta, sigma, gamma, alpha, N, Tend);
disp(Talpha);


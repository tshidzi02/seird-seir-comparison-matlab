function [rmse_seir, rmse_seird, aic_seir, aic_seird] = graphs( ...
    obs, t_seir, y_seir, t_seird, y_seird, sigma)
% GRAPHs  Plot SEIR vs SEIRD incidence and compute metrics.
% Usage:
%   [rmseS, rmseD, aicS, aicD] = graphs(obs, t_seir, y_seir, t_seird, y_seird, sigma)

    % ---------- Basic input validation ----------
    if nargin < 6
        error('graphs:NotEnoughInputs', ...
            'Call as graphs(obs, t_seir, y_seir, t_seird, y_seird, sigma).');
    end
    validateattributes(obs,    {'numeric'},{'vector','nonnegative'});
    validateattributes(t_seir,{'numeric'},{'vector','increasing','nonempty'});
    validateattributes(t_seird,{'numeric'},{'vector','increasing','nonempty'});
    if size(y_seir,2) < 2,  error('y_seir must have at least 2 columns [S E I R].');  end
    if size(y_seird,2) < 2, error('y_seird must have at least 2 columns [S E I R D].'); end
    validateattributes(sigma, {'numeric'},{'scalar','positive'});

    % ---------- Build daily grid from model times (robust to row/col) ----------
    t_seir  = t_seir(:);    % column
    t_seird = t_seird(:);   % column
    T = floor(max([t_seir(end); t_seird(end)]));   % <-- semicolon inside max list avoids size issues
    t_daily = (0:T)';                              

    % ---------- Predicted incidence = sigma * E(t) ----------
    E_seir  = interp1(t_seir,  y_seir(:,2),  t_daily, 'pchip', 'extrap');
    E_seird = interp1(t_seird, y_seird(:,2), t_daily, 'pchip', 'extrap');
    inc_seir  = sigma * E_seir;
    inc_seird = sigma * E_seird;

    % ---------- Align observed with daily grid ----------
    obs      = obs(:);
    L = min([numel(obs), numel(t_daily), numel(inc_seir), numel(inc_seird)]);
    obs      = obs(1:L);
    t_daily  = t_daily(1:L);
    inc_seir = inc_seir(1:L);
    inc_seird= inc_seird(1:L);

    % ---------- Metrics ----------
    n = L;
    k_seir  = 3;    % beta, sigma, gamma
    k_seird = 4;    % beta, sigma, gamma, alpha

    rmse_seir  = sqrt(mean((obs - inc_seir ).^2));
    rmse_seird = sqrt(mean((obs - inc_seird).^2));
    sse_seir   = sum((obs - inc_seir ).^2);
    sse_seird  = sum((obs - inc_seird).^2);
    aic_seir   = 2*k_seir  + n*log(sse_seir /n);
    aic_seird  = 2*k_seird + n*log(sse_seird/n);

    % ---------- Figure 1: Incidence overlay ----------
    figure('Color','w'); hold on;
    plot(t_daily, obs,      'k.', 'DisplayName','Observed', 'MarkerSize', 18);
    plot(t_daily, inc_seir, 'b-' , 'LineWidth',5, 'DisplayName','SEIR');
    plot(t_daily, inc_seird,'r--', 'LineWidth',5, 'DisplayName','SEIRD');
    xlabel('Day'); ylabel('Daily incidence (cases)');
    title(sprintf('Incidence (RMSE SEIR=%.1f, SEIRD=%.1f)', rmse_seir, rmse_seird));
    legend('Location','northwest'); grid on;

    % ---------- Figure 2: Residuals ----------
    res_seir  = obs - inc_seir;
    res_seird = obs - inc_seird;

    figure('Color','w'); tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; plot(t_daily,res_seir,'b-', 'LineWidth',5); yline(0,'k:');
    title('SEIR residuals'); xlabel('Day'); ylabel('Residual'); grid on;
    nexttile; plot(t_daily,res_seird,'r-', 'LineWidth',5); yline(0,'k:');
    title('SEIRD residuals'); xlabel('Day'); ylabel('Residual'); grid on;
end

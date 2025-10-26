function T = sensitivity_alpha(obs, beta, sigma, gamma, alpha, N, Tend, y0)
    if nargin < 8 || isempty(y0), y0 = [N-10;5;5;0;0]; end
    SEIRD = @(a) @(t,y)[-beta*y(1)*y(3)/N; beta*y(1)*y(3)/N - sigma*y(2); ...
                        sigma*y(2) - (gamma+a)*y(3); gamma*y(3); a*y(3)];
    alphas = sort(unique([0 0.001 0.002 0.003 0.01 alpha]));
    t_daily = (0:Tend)'; colors = lines(numel(alphas));
    peak_day = zeros(size(alphas)); rmse = zeros(size(alphas));
    figure('Color','w'); hold on;
    for i=1:numel(alphas)
        [t,y] = ode45(SEIRD(alphas(i)),[0 Tend],y0,odeset('RelTol',1e-7,'AbsTol',1e-9));
        E = interp1(t,y(:,2),t_daily,'pchip'); inc = sigma*max(E,0);
        L = min(numel(inc),numel(obs)); r = obs(1:L)-inc(1:L);
        rmse(i) = sqrt(mean(r.^2)); [~,p]=max(inc); peak_day(i)=t_daily(p);
        plot(t_daily,inc,'LineWidth',3,'Color',colors(i,:),'DisplayName',sprintf('\\alpha=%.3f',alphas(i)));
        % after your plot(...) lines for the Sensitivity-to-α figure:

        set(gcf,'Units','pixels','Position',[100 100 1200 900]);  % wide, tall
        ax = gca;                        % the axes used by this figure
    ax.XTick = 0:10:180;             % tick locations on x
    ax.YTick = 0:2000:25000;         % tick locations on y
    % --- make the window big
    
    
    set(ax,'Position',[0.12 0.12 0.84 0.82]);   % [left bottom width height]
    
    % --- clearer lines & labels
    set(findall(ax,'Type','Line'),'LineWidth',3);   % thicken all curves
    ax.FontSize = 16;
    ax.LineWidth = 1.5;
    legend('Location','northwest','Box','on');     % keep legend off the curves
    
    % --- lighter grid so curves pop
    grid on; grid minor
    ax.GridAlpha = 0.25; ax.MinorGridAlpha = 0.15;
    
    % --- sensible ticks (adjust to taste)
    ax.XTick = 0:10:180;
    ax.YAxis.Exponent = 0; ytickformat('%,.0f');
    
    % --- optional: zoom around the peak to separate curves
    xlim([85 130]);           % focus on the peak window
    ylim([0 24000]);          % or tighten with axis tight/padded
    
    % optional: avoid 10^4 scaling and format nicely
    ax.YAxis.Exponent = 0;
    ytickformat('%,.0f');

    end
    xlabel('Day'); ylabel('Daily incidence'); title('Sensitivity to \alpha (SEIRD)'); legend('Location','northwest'); grid on;
    figure('Color','w');
    subplot(1,2,1); 
    plot(alphas,peak_day,'o-','LineWidth',3); xlabel('\alpha'); ylabel('Peak day'); grid on;
    ax = gca;
    L  = findobj(ax,'Type','Line');           % all line objects in the axes
    for k = 1:numel(L)
        lw = max(2, L(k).LineWidth);          % ensure a decent baseline
        L(k).LineWidth = lw;                  % keep/force your line thickness
        L(k).Marker = 'o';                    % circle markers
        L(k).MarkerSize = 3.2*lw;             % <-- proportional size (try 3–4×)
        L(k).MarkerEdgeColor = L(k).Color;    % match edge to line color
        L(k).MarkerFaceColor = 'w';           % or 'auto' to fill with the line color
    end

    subplot(1,2,2); 
    plot(alphas,rmse,'o-','LineWidth',3); xlabel('\alpha'); ylabel('RMSE'); grid on;
    ax = gca;
    L  = findobj(ax,'Type','Line');           % all line objects in the axes
    for k = 1:numel(L)
        lw = max(2, L(k).LineWidth);          % ensure a decent baseline
        L(k).LineWidth = lw;                  % keep/force your line thickness
        L(k).Marker = 'o';                    % circle markers
        L(k).MarkerSize = 3.2*lw;             % <-- proportional size (try 3–4×)
        L(k).MarkerEdgeColor = L(k).Color;    % match edge to line color
        L(k).MarkerFaceColor = 'w';           % or 'auto' to fill with the line color
    end

    T = table(alphas(:),peak_day(:),rmse(:),'VariableNames',{'alpha','peak_day','RMSE'});
    
end

function plot_cumulative(obs, t_seir, y_seir, t_seird, y_seird, sigma, Tend)
    t_daily = (0:Tend)'; 
    E1 = interp1(t_seir,  y_seir(:,2),  t_daily, 'pchip');
    E2 = interp1(t_seird, y_seird(:,2), t_daily, 'pchip');
    inc1 = sigma*max(E1,0); inc2 = sigma*max(E2,0);
    L = min([numel(obs), numel(t_daily)]);
    cum_obs = cumsum(obs(1:L)); cum1 = cumsum(inc1(1:L)); cum2 = cumsum(inc2(1:L));
   figure('Color','w'); hold on;

% Observed: make the dots bigger (LineWidth has no effect on '.')
hObs  = plot(t_daily(1:L), cum_obs, 'k.', ...
             'MarkerSize', 18, 'DisplayName','Observed');

% Model curves: thicken the lines
hSEIR = plot(t_daily(1:L), cum1, 'b-', ...
             'LineWidth', 4, 'DisplayName','SEIR');
hSEIRD= plot(t_daily(1:L), cum2, 'r--', ...
             'LineWidth', 4, 'DisplayName','SEIRD');
 legend('Observed','SEIR','SEIRD','Location','northwest'); grid on;
  xlabel('Day'); ylabel('Cumulative cases'); title('Cumulative incidence comparison');
set(gca,'FontSize',14,'LineWidth',1.2);

  
end


ax = gca;

% Make all non-dot lines thick
h = findobj(ax,'Type','Line');
for k = 1:numel(h)
    mk = get(h(k),'Marker');
    if isequal(mk,'.') || isequal(mk,'none')
        % dot-only or no marker: enlarge dots if present
        if isequal(mk,'.'), set(h(k),'MarkerSize',18); end
    else
        set(h(k),'LineWidth',4);    % thicken lines
    end
end

% Crisp rendering and readable axes
set(gcf,'GraphicsSmoothing','off');   % avoid anti-aliased “thin” look
set(ax,'FontSize',14,'LineWidth',1.2);
legend('Location','northwest'); grid on


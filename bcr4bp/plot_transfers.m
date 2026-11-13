function figs = plot_transfers(res, halo, prm, fig_dir, picks)
% Plot sample transfers (Earth perigee to halo) in the Earth-Moon and Sun-B1 rotating frames, and their Earth-Moon Jacobi constant over time
% Inputs: res - struct from grid_search; halo - halo struct; prm - struct from bcr4bp_params; fig_dir - folder for the PDFs (optional); picks - rows [Sun angle (deg), halo phase tau, branch p] of viable grid points, plotted as T1, T2, ... (optional, default shortest and median time of flight per branch)
% Outputs: figs - figure handles [transfers, Jacobi constant along them]
    if nargin < 4, fig_dir = ''; end
    if nargin < 5, picks = default_picks(res); end
    mu = prm.mu;  L = prm.lstar;
    figs = figure('Units', 'inches', 'Position', [1 1 14 6.5]);  col = plot_style(figs);
    tl  = tiledlayout(figs, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    a1  = nexttile(tl);  hold(a1, 'on');  axis(a1, 'equal');
    a2  = nexttile(tl);  hold(a2, 'on');  axis(a2, 'equal');
    c   = [col.blue; col.vermillion; col.green; col.orange; col.purple; col.skyblue];
    opts = odeset('RelTol', res.tol, 'AbsTol', res.tol);
    hl  = gobjects(size(picks, 1), 1);
    [td, CJ, dE, sl] = deal(cell(size(picks, 1), 1));  % per transfer: days since departure, C_EM, distance from Earth, label
    rmax = 0;

    for n = 1:size(picks, 1)
        i = find(res.th_deg == picks(n, 1));  j = find(abs(res.taus - picks(n, 2)) < 1e-9);  b = (3 - picks(n, 3))/2;
        if isempty(i) || isempty(j) || ~res.viable(i, j, b), error('Pick %d is not a viable grid point', n); end
        th0 = res.th_deg(i)*pi/180;
        x0  = arrival_states(halo, res.taus(j), picks(n, 3), prm);
        tp  = -res.tof_d(i, j, b)*86400/prm.tstar;
        [t, Y] = ode113(@(t, X) eom_bcr4bp6(t, X, prm, th0), linspace(0, tp, 3000), x0, opts);
        rE = (Y(:, 1:3) - [-mu 0 0])*L/1e5;                         % Earth-centered, 10^5 km
        rS = em_to_sun_frame(Y(:, 1:3), th0 + prm.thdot*t)*L/1e6;  % B1-centered, 10^6 km
        rmax = max(rmax, max(abs(rS(:, 1:2)), [], 'all'));
        hl(n) = plot(a1, rE(:, 1), rE(:, 2), 'Color', c(n, :), 'LineWidth', 1.8);
        plot(a2, rS(:, 1), rS(:, 2), 'Color', c(n, :), 'LineWidth', 1.8);
        plot(a1, rE(end, 1), rE(end, 2), 'o', 'MarkerFaceColor', c(n, :), 'MarkerEdgeColor', 'k');
        plot(a2, rS(end, 1), rS(end, 2), 'o', 'MarkerFaceColor', c(n, :), 'MarkerEdgeColor', 'k');
        set(hl(n), 'DisplayName', sprintf('T%d', n));
        td{n} = (t - tp)*prm.tstar/86400;
        CJ{n} = jacobi_EM(Y, prm);
        dE{n} = vecnorm(rE, 2, 2)/10;                                % 10^6 km
        sl{n} = sprintf('T%d', n);
    end

    ang = linspace(0, 2*pi, 300);
    xH  = halo.sol.y(1:3, :).';
    plot(a1, L/1e5*cos(ang), L/1e5*sin(ang), '--', 'Color', [.6 .6 .6], 'LineWidth', 1);  % Moon's orbit
    plot(a1, (xH(:, 1) + mu)*L/1e5, xH(:, 2)*L/1e5, 'Color', col.black, 'LineWidth', 2);
    hE = plot(a1, 0, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', col.skyblue, 'MarkerEdgeColor', 'k', 'DisplayName', 'Earth');
    hM = plot(a1, L/1e5, 0, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'k', 'DisplayName', 'Moon');
    title(a1, 'Earth-Moon rotating frame (Earth-centered)');
    xlabel(a1, '$x$ ($10^5$ km)');  ylabel(a1, '$y$ ($10^5$ km)');

    plot(a2, (1 - mu)*L/1e6*cos(ang), (1 - mu)*L/1e6*sin(ang), '--', 'Color', [.6 .6 .6], 'LineWidth', 1);
    plot(a2, 0, 0, 'o', 'MarkerSize', 7, 'MarkerFaceColor', col.skyblue, 'MarkerEdgeColor', 'k');
    xline(a2, 0, 'Color', [.75 .75 .75]);  yline(a2, 0, 'Color', [.75 .75 .75]);
    xl = 1.1*rmax;  xlim(a2, [-xl xl]);  ylim(a2, [-xl xl]);
    title(a2, 'Sun-$B_1$ rotating frame (Sun toward $-\hat{x}$)');
    xlabel(a2, '$x$ ($10^6$ km)');  ylabel(a2, '$y$ ($10^6$ km)');
    lg = legend([hl; hE; hM], 'Orientation', 'horizontal');  lg.Layout.Tile = 'south';

    % Earth-Moon Jacobi constant along each transfer (Parker & Anderson Fig. 3-45)
    figs(2) = figure('Units', 'inches', 'Position', [1 1 9 7.5]);  plot_style(figs(2));
    tiledlayout(figs(2), 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    b1 = nexttile;  hold(b1, 'on');  grid(b1, 'on');
    b2 = nexttile;  hold(b2, 'on');  grid(b2, 'on');
    for n = 1:numel(td)
        plot(b1, td{n}, CJ{n}, 'Color', c(n, :), 'DisplayName', sl{n});
        plot(b2, td{n}, dE{n}, 'Color', c(n, :), 'HandleVisibility', 'off');
    end
    Ch = jacobi_EM(halo.sol.y(1:6, 1).', prm);
    yline(b1, Ch, '--', 'Halo', 'Color', col.black, 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'FontSize', 15, 'HandleVisibility', 'off');
    ylim(b1, [0 1.15*Ch]);  ylabel(b1, '$C_{EM}$');  legend(b1, 'Location', 'southeast');
    ylabel(b2, 'Distance from Earth ($10^6$ km)');  xlabel(b2, 'Time since departure (days)');
    linkaxes([b1 b2], 'x');  xlim(b1, [0 max(cellfun(@max, td))]);

    names = {'grid_sample_transfers', 'grid_transfer_jacobi'};
    if ~isempty(fig_dir)
        if ~isfolder(fig_dir), mkdir(fig_dir); end
        for k = 1:2, exportgraphics(figs(k), fullfile(fig_dir, [names{k} '.pdf']), 'ContentType', 'vector'); end
    end
end

function picks = default_picks(res)
% Shortest and median time-of-flight viable points of each branch
% Inputs: res - struct from grid_search
% Outputs: picks - rows [Sun angle (deg), halo phase tau, branch p]
    picks = zeros(0, 3);  th = res.th_deg(:);  ta = res.taus(:);
    for b = 1:2
        T = res.tof_d(:, :, b);  T(~res.viable(:, :, b)) = NaN;
        if all(isnan(T(:))), continue; end
        [~, k1] = min(T(:));
        [~, k2] = min(abs(T(:) - median(T(:), 'omitnan')));
        [i, j] = ind2sub(size(T), [k1; k2]);
        picks  = [picks; th(i), ta(j), 3 - 2*[b; b]]; %#ok<AGROW>
    end
end

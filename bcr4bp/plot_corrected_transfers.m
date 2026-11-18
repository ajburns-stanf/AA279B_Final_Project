function fig = plot_corrected_transfers(out, halo, prm, fig_dir)
% Plot a grid-search seed and its corrected two-impulse transfers: the Earth-Moon rotating, Earth-centered inertial and Sun-B1 rotating frames (top row),
% and zooms on the LEO departure (Earth-centered inertial) and the halo arrival (Earth-Moon rotating) (bottom row)
% Inputs: out - struct from correction_main (seed, z0, seed_info, newton, fmincon); halo - halo struct; prm - struct from bcr4bp_params; fig_dir - folder for the PNG (optional, no export if empty)
% Outputs: fig - figure handle
    if nargin < 4, fig_dir = ''; end
    mu = prm.mu;  L = prm.lstar;  day = 86400/prm.tstar;
    asp = 1.75;  % width/height of the zoom windows, to fill their tiles
    fig = figure('Units', 'inches', 'Position', [0.5 0.5 20 13]);  col = plot_style(fig);
    tl  = tiledlayout(fig, 2, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
    aR  = nexttile(tl, 1, [1 2]);  % Earth-Moon rotating, Earth-centered
    aI  = nexttile(tl, 3, [1 2]);  % Earth-centered inertial
    aS  = nexttile(tl, 5, [1 2]);  % Sun-B1 rotating
    aD  = nexttile(tl, 7, [1 3]);  % departure zoom, Earth-centered inertial
    aA  = nexttile(tl, 10, [1 3]); % arrival zoom, Earth-Moon rotating, Moon-centered
    ax  = [aR aI aS aD aA];
    hold(ax, 'on');  axis(aR, 'equal');  axis(aI, 'equal');  axis(aS, 'equal');
    opts = odeset('RelTol', prm.dc.tol, 'AbsTol', prm.dc.tol);

    % Seed (grid arc, no HOI), Newton (gap closed), fmincon (minimum TLI + HOI)
    I  = {out.seed_info, out.newton.info, out.fmincon.info};
    Z  = {out.z0, out.newton.z, out.fmincon.z};
    nm = {'Seed', 'Newton', 'fmincon'};
    ls = {'--', '-', '-'};
    c  = [col.black; col.blue; col.vermillion];
    h  = gobjects(0);  dep = zeros(0, 2);
    for n = find([true, true, out.fmincon.converged])
        z  = Z{n};
        xf = arrival_states(halo, z(5), out.seed.p, prm) - [0; 0; 0; z(1:3)];
        t  = sort(unique([linspace(0, z(4), 4000), linspace(z(4) + day, z(4), 2000), linspace(z(4) + day/8, z(4), 1000)]), 'descend');  % dense after departure
        [t, Y] = ode113(@(t, X) eom_bcr4bp6(t, X, prm, z(6)), t, xf, opts);

        rE = (Y(:, 1:3) - [-mu 0 0])*L;                                       % Earth-centered rotating, km
        rM = (Y(:, 1:3) - [1 - mu 0 0])*L;                                    % Moon-centered rotating, km
        a  = t - z(4);                                                        % frame angle since departure (rad, frame rate 1)
        rI = [cos(a).*rE(:, 1) - sin(a).*rE(:, 2), sin(a).*rE(:, 1) + cos(a).*rE(:, 2), rE(:, 3)];  % Earth-centered inertial, km
        rS = em_to_sun_frame(Y(:, 1:3), z(6) + prm.thdot*t)*L;               % B1-centered Sun-B1 rotating, km

        lbl = sprintf('%s: %.4f km/s', nm{n}, I{n}.dv_kms);
        if n == 1, lbl = sprintf('Seed: $h_p$ = %.1f km', I{n}.hp_km); end
        sty = {ls{n}, 'Color', c(n, :), 'LineWidth', 1.6};
        mA  = {'o', 'MarkerSize', 7, 'MarkerFaceColor', c(n, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off'};  % halo arrival (the arc runs backward from it)
        mD  = {'s', 'MarkerSize', 8, 'MarkerFaceColor', c(n, :), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off'};  % departure
        h(end+1) = plot(aR, rE(:, 1)/1e5, rE(:, 2)/1e5, sty{:}, 'DisplayName', lbl); %#ok<AGROW>
        plot(aI, rI(:, 1)/1e6, rI(:, 2)/1e6, sty{:});
        plot(aS, rS(:, 1)/1e6, rS(:, 2)/1e6, sty{:});
        plot(aD, rI(:, 1)/1e3, rI(:, 2)/1e3, sty{:});
        plot(aA, rM(:, 1)/1e3, rM(:, 2)/1e3, sty{:});
        k   = a < day & vecnorm(rI, 2, 2) < 4.5*prm.rLEO*L;  % departure leg near the Earth, frames the zoom
        dep = [dep; rI(k, 1:2)/1e3]; %#ok<AGROW>
        plot(aR, rE(1, 1)/1e5, rE(1, 2)/1e5, mA{:});
        plot(aI, rI(1, 1)/1e6, rI(1, 2)/1e6, mA{:});
        plot(aS, rS(1, 1)/1e6, rS(1, 2)/1e6, mA{:});
        plot(aA, rM(1, 1)/1e3, rM(1, 2)/1e3, mA{:});
        plot(aD, rI(end, 1)/1e3, rI(end, 2)/1e3, mD{:});
    end
    h(end+1) = plot(aR, NaN, NaN, 's', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Departure');
    h(end+1) = plot(aR, NaN, NaN, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Halo arrival');

    ang  = linspace(0, 2*pi, 300);
    xH   = halo.sol.y(1:3, :).';
    body = @(ax, x, y, ms, fc) plot(ax, x, y, 'o', 'MarkerSize', ms, 'MarkerFaceColor', fc, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    grey = {':', 'Color', [.6 .6 .6], 'LineWidth', 1, 'HandleVisibility', 'off'};
    haloline = {'Color', col.green, 'LineWidth', 2, 'HandleVisibility', 'off'};

    % Earth-Moon rotating
    plot(aR, L/1e5*cos(ang), L/1e5*sin(ang), grey{:});
    plot(aR, (xH(:, 1) + mu)*L/1e5, xH(:, 2)*L/1e5, haloline{:});
    body(aR, 0, 0, 9, col.skyblue);  body(aR, L/1e5, 0, 5, [.5 .5 .5]);
    title(aR, 'Earth-Moon rotating (Earth-centered)');
    xlabel(aR, '$x$ ($10^5$ km)');  ylabel(aR, '$y$ ($10^5$ km)');

    % Earth-centered inertial
    plot(aI, L/1e6*cos(ang), L/1e6*sin(ang), grey{:});
    body(aI, 0, 0, 7, col.skyblue);
    title(aI, 'Earth-centered inertial ($\hat{x}$: Earth-Moon line at departure)');
    xlabel(aI, '$X$ ($10^6$ km)');  ylabel(aI, '$Y$ ($10^6$ km)');

    % Sun-B1 rotating
    plot(aS, (1 - mu)*L/1e6*cos(ang), (1 - mu)*L/1e6*sin(ang), grey{:});
    xline(aS, 0, 'Color', [.75 .75 .75], 'HandleVisibility', 'off');  yline(aS, 0, 'Color', [.75 .75 .75], 'HandleVisibility', 'off');
    body(aS, 0, 0, 7, col.skyblue);
    xl = xlim(aS);  text(aS, xl(1), 0, '$\leftarrow$ Sun', 'VerticalAlignment', 'bottom', 'FontSize', 15);
    title(aS, 'Sun-$B_1$ rotating (Sun toward $-\hat{x}$)');
    xlabel(aS, '$x$ ($10^6$ km)');  ylabel(aS, '$y$ ($10^6$ km)');

    % Departure zoom: Earth disk and the LEO radius
    RE = prm.RE*L/1e3;  rL = prm.rLEO*L/1e3;
    fill(aD, RE*cos(ang), RE*sin(ang), col.skyblue, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(aD, rL*cos(ang), rL*sin(ang), '-.', 'Color', [.35 .35 .35], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(aD, rL*cos(-3*pi/4), rL*sin(-3*pi/4), sprintf('LEO %d km  ', prm.hLEO_km), 'Color', [.35 .35 .35], 'FontSize', 15, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    bx = [min([dep(:, 1); -rL]), max([dep(:, 1); rL])];  by = [min([dep(:, 2); -rL]), max([dep(:, 2); rL])];
    w  = 1.08*max(diff(bx), asp*diff(by));  % box around the Earth and the departure legs, widened to the zoom aspect
    set(aD, 'XLim', mean(bx) + w/2*[-1 1], 'YLim', mean(by) + w/asp/2*[-1 1], 'DataAspectRatio', [1 1 1]);
    title(aD, 'Departure (Earth-centered inertial)');
    xlabel(aD, '$X$ ($10^3$ km)');  ylabel(aD, '$Y$ ($10^3$ km)');

    % Arrival zoom: halo and the Moon
    xHM = (xH(:, 1) - (1 - mu))*L/1e3;  yHM = xH(:, 2)*L/1e3;
    plot(aA, xHM, yHM, haloline{:});
    body(aA, 0, 0, 9, [.5 .5 .5]);
    text(aA, 0, 0, 'Moon  ', 'HorizontalAlignment', 'right', 'FontSize', 15);
    hy = 1.2*max(abs(yHM));  w = 2*asp*hy;  % Moon at the left edge, the approach from beyond the halo fills the rest
    set(aA, 'XLim', -0.08*w + [0 w], 'YLim', hy*[-1 1], 'DataAspectRatio', [1 1 1]);
    title(aA, 'Halo arrival (Earth-Moon rotating, Moon-centered)');
    xlabel(aA, '$x - x_M$ ($10^3$ km)');  ylabel(aA, '$y$ ($10^3$ km)');

    lg = legend(aR, h, 'Orientation', 'horizontal');  lg.Layout.Tile = 'south';

    if ~isempty(fig_dir)
        if ~isfolder(fig_dir), mkdir(fig_dir); end
        exportgraphics(fig, fullfile(fig_dir, sprintf('corrected_th%g_tau%.2f_p%+d.png', out.seed.thS_deg, out.seed.tau, out.seed.p)), 'Resolution', 200);
    end
end

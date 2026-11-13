%% BCR4BP checks that produce figures
% Each section makes one figure and prints the number it checks; run after AA279B_FinalProject.m (uses its halo and prm)
if ~exist('halo', 'var') || ~exist('prm', 'var'), error('Run AA279B_FinalProject.m first: this script uses its halo and prm'); end
here = fileparts(mfilename('fullpath'));
if isempty(here), here = fileparts(which('check_bcr4bp_figures')); end  % fallback if mfilename is empty
addpath(here, fileparts(here));  % bcr4bp/ and repo root (plot_style)
save_figs = true;                % false skips PDF export
fig_dir   = fullfile(fileparts(fileparts(here)), 'figures');  % next to the repo, never inside it (report figures)

p0   = prm;  p0.muS = 0;         % muS = 0 is the Earth-Moon CR3BP
mu   = prm.mu;  L = prm.lstar;  day = 86400/prm.tstar;
tol  = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
opts = odeset(tol, 'Events', @(t, X) hit_body(t, X, prm));

% Corrected halo from the main script: CR3BP solution over one period
solW = halo.sol;  P = halo.P;
X0   = solW.y(1:6, 1);
xH   = solW.y(1:3, :).';         % halo positions over one period


%% CR3BP halo IC propagated in the BCR4BP (Scheuerle Fig. 3.10 p. 102)
th5 = [0 90 180 270];
ls5 = {'-', '-', '--', '--'};                % 180 and 270 nearly overlay 0 and 90 (tidal symmetry)
tt  = linspace(0, 4*P, 4000);
rH  = @(t) deval(solW, mod(t, P), 1:3).';    % periodic CR3BP halo
mx  = @(x) (x - (1 - mu))*L/1e3;             % Moon-centered, 10^3 km
my  = @(y) y*L/1e3;

fig = figure('Units', 'inches', 'Position', [1 1 14 6]);  col = plot_style(fig);
c5  = [col.blue; col.vermillion; col.green; col.orange];
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile; hold on; axis equal; grid on;
ax2 = nexttile; hold on; set(ax2, 'YScale', 'log'); grid on;

[t, Y] = ode113(@(t, X) eom_bcr4bp6(t, X, p0, 0), tt, X0, opts);
dev = vecnorm(Y(:, 1:3) - rH(t), 2, 2)*L;
fprintf('Halo drift: CR3BP         deviation passes 100 km at %.1f d\n', first_time(t, dev, 100)/day);
for k = 1:numel(th5)
    [t, Y] = ode113(@(t, X) eom_bcr4bp6(t, X, prm, th5(k)*pi/180), tt, X0, opts);
    dev = vecnorm(Y(:, 1:3) - rH(t), 2, 2)*L;
    k1  = t <= P;                            % first revolution only
    plot(ax1, mx(Y(k1, 1)), my(Y(k1, 2)), ls5{k}, 'Color', c5(k, :));
    plot(ax2, t/day, dev, ls5{k}, 'Color', c5(k, :), 'DisplayName', sprintf('$\\theta_{s,0} = %d^\\circ$', th5(k)));
    fprintf('Halo drift: thS0 = %3d deg, deviation passes 100 km at %.1f d and 1000 km at %.1f d\n', ...
        th5(k), first_time(t, dev, 100)/day, first_time(t, dev, 1000)/day);
end
hH = plot(ax1, mx(xH(:, 1)), my(xH(:, 2)), 'Color', col.black, 'LineWidth', 3, 'DisplayName', 'CR3BP halo');
plot(ax1, 0, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'k');
title(ax1, 'First revolution in Moon-centered frame');
xlabel(ax1, '$x$ ($10^3$ km)');  ylabel(ax1, '$y$ ($10^3$ km)');
yline(ax2, 100, '-.', 'Color', [.5 .5 .5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
title(ax2, 'Distance from halo over time');
xlabel(ax2, '$t$ (days)');  ylabel(ax2, 'Distance from halo (km)');
xlim(ax2, [0 4*P/day]);  ylim(ax2, [1 3e6]);
lg = legend([hH; ax2.Children(end:-1:1)], 'Location', 'eastoutside');  lg.Layout.Tile = 'east';  % one legend for both panels
save_pdf(fig, 'halo_in_bcr4bp', save_figs, fig_dir);


%% Stable manifold in the BCR4BP, Sun-B1 frame (Parker & Anderson Fig. 3-35 p. 169, Scheuerle Fig. 4.53 p. 175)
% Grid-search arrival states (arrival_states): 100 halo phases, both branches, 100 km shift, 180 d backward
taus = (0:99)/100;
Xs   = {arrival_states(halo, taus, 1, prm), arrival_states(halo, taus, -1, prm)};
th6  = [0 45 90 135];
tt6  = linspace(0, -180*day, 3000);
cm   = (1 - mu)*L/1e6*[cos(linspace(0, 2*pi, 200)); sin(linspace(0, 2*pi, 200))];

fig = figure('Units', 'inches', 'Position', [1 1 10 11.5]);  col = plot_style(fig);
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(th6)
    th0 = th6(k)*pi/180;
    ax  = nexttile; hold on; axis equal;
    viab = {};  nmoon = 0;
    for j = 1:numel(taus)
        for p = [1 -1]
            [t, Y, ~, ~, ie] = ode113(@(t, X) eom_bcr4bp6(t, X, prm, th0), tt6, Xs{(3 - p)/2}(:, j), opts);
            rS = em_to_sun_frame(Y(:, 1:3), th0 + prm.thdot*t)*L/1e6;
            plot(ax, rS(:, 1), rS(:, 2), 'Color', [c_branch(col, p) 0.25], 'LineWidth', 0.5);
            [rp, kp] = min(vecnorm(Y(:, 1:3) - [-mu 0 0], 2, 2));
            hp = (rp - prm.RE)*L;
            if hp <= prm.hLEO_km, viab{end+1} = rS(1:kp, :); end %#ok<SAGROW> % viable as in grid_search (Earth impact ends at hp = 0): halo back to perigee
            nmoon = nmoon + any(ie == 2);
        end
    end
    for j = 1:numel(viab), plot(ax, viab{j}(:, 1), viab{j}(:, 2), 'Color', col.black, 'LineWidth', 1.2); end
    hO = plot(ax, cm(1, :), cm(2, :), '--', 'Color', [.3 .3 .3], 'LineWidth', 1.8, 'DisplayName', 'Moon''s orbit');
    hE = plot(ax, 0, 0, 'o', 'MarkerSize', 7, 'MarkerFaceColor', col.skyblue, 'MarkerEdgeColor', 'k', 'DisplayName', 'Earth');
    xline(ax, 0, 'Color', [.7 .7 .7]);  yline(ax, 0, 'Color', [.7 .7 .7]);
    xlim(ax, [-2 2]);  ylim(ax, [-2 2]);
    title(ax, sprintf('$\\theta_{s,0} = %d^\\circ$', th6(k)));
    xlabel(ax, '$x$ ($10^6$ km)');  ylabel(ax, '$y$ ($10^6$ km)');
    fprintf('Manifold: thS0 = %3d deg, viable arcs (perigee alt <= %d km or Earth impact): %d / %d, Moon impacts: %d\n', ...
        th6(k), prm.hLEO_km, numel(viab), 2*numel(taus), nmoon);
end
title(tl, 'Sun-$B_1$ rotating frame (Sun toward $-\hat{x}$)', 'Interpreter', 'latex', 'FontSize', 20);
hb = [plot(ax, NaN, NaN, 'Color', col.blue, 'LineWidth', 2.5, 'DisplayName', 'Exterior');     % legend-only copies, bold and opaque
      plot(ax, NaN, NaN, 'Color', col.vermillion, 'LineWidth', 2.5, 'DisplayName', 'Interior');
      plot(ax, NaN, NaN, 'Color', col.black, 'LineWidth', 2.5, 'DisplayName', sprintf('Perigee alt. $\\leq %d$ km', prm.hLEO_km)); hO; hE];
lg = legend(hb, 'Orientation', 'horizontal');
lg.Layout.Tile = 'south';
save_pdf(fig, 'manifold_bcr4bp', save_figs, fig_dir, 'image');


function [g, term, dirn] = hit_body(~, X, prm)
% Event: trajectory reaches the Earth or Moon surface
% Inputs: X - state (position in entries 1-3), prm - struct from bcr4bp_params
% Outputs: g - distances above the surfaces, term - stop flags, dirn - crossing directions
    g    = [norm(X(1:3) - [-prm.mu; 0; 0]) - prm.RE; norm(X(1:3) - [1 - prm.mu; 0; 0]) - prm.RM];
    term = [1; 1];  dirn = [0; 0];
end

function t1 = first_time(t, d, level)
% First time a distance history exceeds a level
% Inputs: t - times, d - distances, level - threshold
% Outputs: t1 - first time with d > level (NaN if never)
    k  = find(d > level, 1);
    t1 = NaN;
    if ~isempty(k), t1 = t(k); end
end

function c = c_branch(col, p)
% Color for a manifold branch
% Inputs: col - struct from plot_style, p - branch sign (+1 or -1)
% Outputs: c - RGB color
    if p > 0, c = col.blue; else, c = col.vermillion; end
end

function save_pdf(fig, name, save_figs, fig_dir, ctype)
% Export a figure to PDF when saving is on
% Inputs: fig - figure, name - file name without extension, save_figs - flag, fig_dir - output folder, ctype - 'vector' (default) or 'image' (300 dpi, for dense plots)
% Outputs: none
    if ~save_figs, return; end
    if nargin < 5, ctype = 'vector'; end
    if ~isfolder(fig_dir), mkdir(fig_dir); end
    if strcmp(ctype, 'image')
        exportgraphics(fig, fullfile(fig_dir, [name '.pdf']), 'ContentType', 'image', 'Resolution', 300);
    else
        exportgraphics(fig, fullfile(fig_dir, [name '.pdf']), 'ContentType', 'vector');
    end
end

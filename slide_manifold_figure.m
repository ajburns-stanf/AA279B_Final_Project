%% AA 279B Final Project: CR3BP Stable Manifold Figure for Slides
% Propagates the L2 halo stable manifold in the CR3BP (BCR4BP with the Sun's mass set to zero),
% using the grid-search settings (100 phases, both branches, 100 km shift, 180 days backward),
% and reports the closest approach to Earth. Figure saved outside the repo, like the grid figures.
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'bcr4bp')); % BCR4BP functions

prm  = bcr4bp_params();
prm0 = prm;  prm0.muS = 0;          % CR3BP: Sun's mass removed
col  = plot_style();
mu   = prm.mu;

% -- L2 Halo Orbit + STM (Earth-side IC, already differentially corrected in the main script) --
X0 = [1.0606924537198088; 0; 0.07354804500048147; 0; 0.3422130947784463; 0];
P  = 3.2236501365076409;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
halo.sol = ode113(@(t,w) eom_stm(t, w, prm0), [0 P], [X0; reshape(eye(6), [], 1)], opts);
halo.P   = P;

% -- Stable Eigenvector of the Monodromy Matrix --
M = reshape(halo.sol.y(7:42, end), 6, 6);
[V, Lam] = eig(M);
[~, is]  = min(abs(diag(Lam)));
halo.vs  = real(V(:, is));
fprintf('Periodicity miss %.2e, lambda_u = %.2f\n', norm(halo.sol.y(1:6, end) - X0), max(abs(diag(Lam))));

% -- Manifold Arcs, Grid-Search Settings --
tau  = (0:99)/100;  branches = [+1 -1];   % +1 exterior, -1 interior
tend = -prm.dt_d*86400/prm.tstar;
evo  = odeset('RelTol', prm.tol, 'AbsTol', prm.tol, 'Events', @(t, X) stop_events(t, X, prm0));
arcs = cell(numel(branches), numel(tau));
hmin = inf(numel(branches), numel(tau));  tmin = nan(size(hmin));
for b = 1:numel(branches)
    x0 = arrival_states(halo, tau, branches(b), prm0);
    for k = 1:numel(tau)
        out = propagate_grid_point(x0(:, k), 0, prm0);          % true perigees, same as the grid search
        hmin(b, k) = out.hmin_km;  tmin(b, k) = out.tof_d;
        [~, X] = ode113(@(t, X) eom_bcr4bp6(t, X, prm0, 0), [0 tend], x0(:, k), evo);
        arcs{b, k} = X;
    end
end

[hbest, ib] = min(hmin(:));
[bb, kb]    = ind2sub(size(hmin), ib);
fprintf('Closest CR3BP manifold perigee: %.0f km altitude (tau = %.2f, p = %+d, %.0f days before arrival)\n', ...
        hbest, tau(kb), branches(bb), tmin(bb, kb));
fprintf('Arcs with a perigee below %d km: %d of %d\n', prm.hLEO_km, nnz(hmin <= prm.hLEO_km), numel(hmin));

% -- Plot: Earth-Moon rotating frame, Earth-centered, xy-plane --
s  = prm.lstar/1e5;                       % 10^5 km
ex = @(X) (X(:, 1) + mu)*s;  ey = @(X) X(:, 2)*s;
bcol = {col.blue, col.orange};

fig = figure(1); clf; hold on; grid on; axis equal
set(fig, 'Units', 'inches', 'Position', [1 1 7.5 7.5])
for b = 1:numel(branches)
    for k = 1:numel(tau)
        hB(b) = plot(ex(arcs{b, k}), ey(arcs{b, k}), 'Color', [bcol{b} 0.35], 'LineWidth', 0.6); %#ok<SAGROW>
    end
end
th = linspace(0, 2*pi, 400);
hMo = plot(s*cos(th), s*sin(th), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);      % Moon's orbit
W   = deval(halo.sol, linspace(0, P, 400));
hH  = plot((W(1, :) + mu)*s, W(2, :)*s, 'r', 'LineWidth', 2.5);                        % L2 halo

% Closest approach: circle of that radius about Earth + marker at the perigee
rc  = (hbest + 6378.137)/1e5;
hC  = plot(rc*cos(th), rc*sin(th), ':', 'Color', col.purple, 'LineWidth', 2);
Xb  = arcs{bb, kb};  [~, j] = min(vecnorm(Xb(:, 1:3) - [-mu 0 0], 2, 2));
plot(ex(Xb(j, :)), ey(Xb(j, :)), 'p', 'MarkerSize', 14, 'MarkerFaceColor', col.purple, 'MarkerEdgeColor', 'k')

hE = plot(0, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', col.blue, 'MarkerEdgeColor', 'k');      % Earth
plot(s, 0, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'k')          % Moon

xlim([-6 8]); ylim([-6 6])
xlabel('$x$ ($10^5$ km)'); ylabel('$y$ ($10^5$ km)')
hc_str = regexprep(sprintf('%d', round(hbest, -3)), '(\d)(?=(\d{3})+$)', '$1,');   % e.g. 63,000
legend([hB(1) hB(2) hH hMo hC], {'Exterior', 'Interior', 'L$_2$ halo', 'Moon''s orbit', ...
        ['Closest pass: ' hc_str ' km alt.']}, 'Location', 'southoutside', 'NumColumns', 3)

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');   % outside the repo
if ~isfolder(out_dir), mkdir(out_dir); end
exportgraphics(fig, fullfile(out_dir, 'slide_manifold_cr3bp.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, 'slide_manifold_cr3bp.pdf'), 'ContentType', 'vector');
fprintf('Saved to %s\n', out_dir);


%% ======= Helper Functions =======

% ---- CR3BP State + STM (Phi_dot = A*Phi), A from the pseudo-potential Hessian ----
function wdot = eom_stm(t, w, prm)
    mu = prm.mu;  x = w(1); y = w(2); z = w(3);
    r1 = sqrt((x + mu)^2 + y^2 + z^2);  r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    a = (1 - mu)/r1^3;  b = mu/r2^3;  A1 = 3*(1 - mu)/r1^5;  B1 = 3*mu/r2^5;
    Uxx = 1 - a - b + A1*(x + mu)^2 + B1*(x - 1 + mu)^2;
    Uyy = 1 - a - b + (A1 + B1)*y^2;
    Uzz =   - a - b + (A1 + B1)*z^2;
    Uxy = A1*(x + mu)*y + B1*(x - 1 + mu)*y;
    Uxz = A1*(x + mu)*z + B1*(x - 1 + mu)*z;
    Uyz = (A1 + B1)*y*z;
    A = [zeros(3) eye(3); [Uxx Uxy Uxz; Uxy Uyy Uyz; Uxz Uyz Uzz] [0 2 0; -2 0 0; 0 0 0]];
    Phi  = reshape(w(7:42), 6, 6);
    wdot = [eom_bcr4bp6(t, w(1:6), prm, 0); reshape(A*Phi, [], 1)];
end

% ---- Stop at Earth impact, Moon impact or escape (same limits as the grid search) ----
function [g, term, dirn] = stop_events(~, X, prm)
    g    = [norm(X(1:3) - [-prm.mu; 0; 0]) - prm.RE; norm(X(1:3) - [1 - prm.mu; 0; 0]) - prm.RM; norm(X(1:3)) - prm.r_esc];
    term = [1; 1; 1];
    dirn = [0; 0; 0];
end

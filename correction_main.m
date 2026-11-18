%% AA 279B Final Project: Differential Correction of a Grid-Search Seed
% = Alexander Burns, Kaila Coimbra, Yuji Takubo =
% Two-impulse transfer: tangential TLI from the circular LEO, HOI at the stable-manifold arrival state.
% Loads the grid that AA279B_FinalProject.m saves (results/grid_coarse.mat, next to the repo); without it, propagates only the chosen seed.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, 'bcr4bp'));

prm  = bcr4bp_params();
halo = build_halo(prm);  % the L2 halo and stable eigenvector of Part 1 of AA279B_FinalProject.m

out_dir = fileparts(here);              % folder holding the repo, results saved outside repo
res_dir = fullfile(out_dir, 'results');
fig_dir = fullfile(here, 'figure');     % PNGs in the repo
if ~isfolder(res_dir), mkdir(res_dir); end


%% ========= Part 1 ===============
%   Correct One Grid-Search Seed
% =================================
disp('Part 1: Correct one grid-search seed')

% -- Seed: [Sun angle at arrival (deg), halo phase tau, branch p], a true-perigee grid point (T1-T4 of the main script) --
pick = [10 0.66 1];  % T1; also T2 [180 0.54 -1], T3 [280 0.07 -1]; T4 [25 0.13 -1] has a lunar flyby and is screened out
seed = grid_seed(pick, fullfile(res_dir, 'grid_coarse.mat'), halo, prm);
tag  = sprintf('th%g_tau%.2f_p%+d', pick(1), pick(2), pick(3));

% -- Setup and initial guess: no HOI, the seed's time of flight and epoch --
P  = dc_setup(seed.p, halo, prm);
z0 = [0; 0; 0; -seed.tof_d*86400/prm.tstar; seed.tau; seed.thS_deg*pi/180];

% -- Stage 1: damped Newton closes F = 0 from the seed, moving dv and tp (tau and thS0 stay at the seed) --
tic;
opt_newton = dc_newton(P, z0, P.newton_vars);
if ~opt_newton.converged
    error('correction_main:newton', 'Newton did not converge (|F| = %.2e after %d iterations): pick another seed', ...
        norm(opt_newton.F), opt_newton.iter);
end
fprintf('Stage 1 (Newton):  %.1f s\n', toc);

% -- Stage 2: fmincon minimizes J over the outer variables dv and tau; inside every cost evaluation, --
% -- Newton re-solves the inner variables tp and thS0 so that F = 0 (see dc_optimize)                 --
tic;
opt_fmincon = dc_optimize(P, opt_newton);
fprintf('Stage 2 (fmincon): %.1f s\n\n', toc);

% -- Compare the seed, Newton and fmincon (fmincon's transfer is kept if opt_fmincon.converged) --
seed_eval = P.eval(z0);  % the seed arc itself, for its perigee
out = struct('seed', seed, 'z0', z0, 'seed_info', seed_eval.info, 'newton', opt_newton, 'fmincon', opt_fmincon);
print_summary(out);

% -- Save and plot --
save(fullfile(res_dir, ['corrected_' tag '.mat']), 'out', 'prm');
plot_corrected_transfers(out, halo, prm, fig_dir);



%% ======= Helper Functions =======

% ---- Seed from the grid results, or from one backward propagation if they are missing ----
function seed = grid_seed(pick, res_file, halo, prm)
    seed = struct('thS_deg', pick(1), 'tau', pick(2), 'p', pick(3));
    if isfile(res_file)
        S = load(res_file, 'res');  res = S.res;
        i = find(res.th_deg == pick(1));  j = find(abs(res.taus - pick(2)) < 1e-9);  b = (3 - pick(3))/2;
        if isempty(i) || isempty(j), error('Pick [%g %g %d] is not a grid point', pick); end
        [seed.tof_d, seed.hmin_km, seed.dv_kms] = deal(res.tof_d(i, j, b), res.hmin_km(i, j, b), res.dv_kms(i, j, b));
        fprintf('Seed from %s\n', res_file);
    else
        warning('correction_main:nogrid', '%s not found (run AA279B_FinalProject.m); propagating this seed only', res_file);
        o = propagate_grid_point(arrival_states(halo, pick(2), pick(3), prm), pick(1)*pi/180, prm);
        [seed.tof_d, seed.hmin_km, seed.dv_kms] = deal(o.tof_d, o.hmin_km, o.dv_kms);
    end
    if ~(seed.hmin_km > 0 && isfinite(seed.tof_d))
        error('Pick [%g %g %d] has no true perigee (Earth impact or none in the window)', pick);
    end
    [flyby, rmin_km, tmin_d] = moon_flyby(seed, halo, prm);
    if flyby
        error('Pick [%g %g %d] has a lunar flyby (%.0f km from the Moon, %.1f days before arrival): too sensitive to correct, pick another seed', ...
            pick, rmin_km, tmin_d);
    end
end


% ---- Table: grid estimate, Newton, fmincon ----
function print_summary(out)
    s = out.seed;  N = out.newton;  O = out.fmincon;
    fprintf('Seed: thS = %g deg, tau = %.2f, p = %+d\n', s.thS_deg, s.tau, s.p);
    fprintf('%-8s %10s %10s %10s %12s %11s %9s %9s %9s %6s\n', '', 'TLI km/s', 'HOI m/s', 'Total km/s', 'hp km', 'FPA deg', 'TOF d', 'tau', 'thS deg', 'iter');
    fprintf('%-8s %10.5f %10s %10.5f %12.3f %11s %9.3f %9.5f %9.4f %6s\n', 'Grid', s.dv_kms, '-', s.dv_kms, s.hmin_km, '-', s.tof_d, s.tau, s.thS_deg, '-');
    rows = {'Newton', N; 'fmincon', O};
    for k = 1:2
        I = rows{k, 2}.info;  z = rows{k, 2}.z;
        fprintf('%-8s %10.5f %10.4f %10.5f %12.6f %11.2e %9.3f %9.5f %9.4f %6d\n', rows{k, 1}, I.tli_kms, I.hoi_kms*1e3, I.dv_kms, ...
            I.hp_km, I.fpa_deg, I.tof_d, z(5), z(6)*180/pi, rows{k, 2}.iter);
    end
    fprintf('Grid TLI is the energy patch from the seed perigee (propagate_grid_point); Newton and fmincon are exact at the LEO.\n');
    fprintf('Newton converged: %d (|F| = %.1e); fmincon kept: %d (exitflag %d, %d cost evaluations, %d inner solves failed)\n', ...
        N.converged, norm(N.F), O.converged, O.exitflag, O.evals, O.failed);
end

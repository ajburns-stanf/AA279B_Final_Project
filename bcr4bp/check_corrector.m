%% Differential-corrector checks
% Derivatives against central finite differences, then the two-stage correction (dc_newton, dc_optimize) on the main script's sample transfers T1-T4 with an independent re-propagation
% Self-contained: builds the halo and propagates the seeds itself (no grid file needed)
here = fileparts(mfilename('fullpath'));
if isempty(here), here = fileparts(which('check_corrector')); end
addpath(here, fileparts(here));
prm   = bcr4bp_params();
halo  = build_halo(prm);
picks = [10 0.66 1; 180 0.54 -1; 280 0.07 -1; 25 0.13 -1];  % T1-T4: [Sun angle (deg), tau, p]
seeds = cell(size(picks, 1), 1);
for n = 1:size(picks, 1)
    o = propagate_grid_point(arrival_states(halo, picks(n, 2), picks(n, 3), prm), picks(n, 1)*pi/180, prm);
    seeds{n} = struct('thS_deg', picks(n, 1), 'tau', picks(n, 2), 'p', picks(n, 3), 'tof_d', o.tof_d, 'hmin_km', o.hmin_km, 'dv_kms', o.dv_kms);
end
relerr = @(a, b) norm(a - b)/norm(b);


%% 1. Arrival-state derivative d(xf)/d(tau) (arrival_states second output)
for n = 1:3
    [~, dxf] = arrival_states(halo, seeds{n}.tau, seeds{n}.p, prm);
    h  = 1e-6;
    fd = (arrival_states(halo, seeds{n}.tau + h, seeds{n}.p, prm) - arrival_states(halo, seeds{n}.tau - h, seeds{n}.p, prm))/(2*h);
    fprintf('1. T%d d(xf)/d(tau): relative error %.1e\n', n, relerr(dxf, fd));
end


%% 2. STM and Sun-angle sensitivity of eom_bcr4bp6 over a 10-day backward arc from T1's arrival state
s1 = seeds{1};  th0 = s1.thS_deg*pi/180;  T = -10*86400/prm.tstar;
x0 = arrival_states(halo, s1.tau, s1.p, prm);
o6 = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
fly = @(x, th) deval(ode113(@(t, X) eom_bcr4bp6(t, X, prm, th), [0 T], x, o6), T);
W   = fly([x0; reshape(eye(6), [], 1); zeros(6, 1)], th0);
Phi = reshape(W(7:42), 6, 6);  psi = W(43:48);
fdP = zeros(6);  h = 1e-7;
for k = 1:6
    e = zeros(6, 1);  e(k) = h;
    fdP(:, k) = (fly(x0 + e, th0) - fly(x0 - e, th0))/(2*h);
end
fdS = (fly(x0, th0 + h) - fly(x0, th0 - h))/(2*h);
fprintf('2. STM: relative error %.1e; d(X)/d(thS0): relative error %.1e; 6-state vs 48-state endpoint %.1e\n', ...
    relerr(Phi, fdP), relerr(psi, fdS), norm(fly(x0, th0) - W(1:6)));


%% 3. Residual and cost gradients of correction_obj (full T1-T3 arcs, at each seed)
hz = [1e-10 1e-10 1e-10 1e-8 1e-9 1e-9];  % steps: dv (VU), tp (TU), tau, thS0 (rad)
for n = 1:3
    z = [0; 0; 0; -seeds{n}.tof_d*86400/prm.tstar; seeds{n}.tau; seeds{n}.thS_deg*pi/180];
    [~, dF, ~, dJ] = correction_obj(z, seeds{n}.p, halo, prm);
    [fdF, fdJ] = deal(zeros(2, 6), zeros(1, 6));
    for k = 1:6
        e = zeros(6, 1);  e(k) = hz(k);
        [Fp, ~, Jp] = correction_obj(z + e, seeds{n}.p, halo, prm);
        [Fm, ~, Jm] = correction_obj(z - e, seeds{n}.p, halo, prm);
        fdF(:, k) = (Fp - Fm)/(2*hz(k));  fdJ(k) = (Jp - Jm)/(2*hz(k));
    end
    fprintf('3. T%d dF/dz column errors %s; dJ/dz %.1e\n', n, mat2str(vecnorm(dF - fdF)./vecnorm(fdF), 2), relerr(dJ, fdJ));
end


%% 4. Two-stage correction (dc_newton, then dc_optimize, as in correction_main) on T1-T4, then an independent 6-state re-propagation of each kept transfer
fprintf('\n4. %-3s %9s %9s %9s %11s | %10s %8s %10s | %10s %9s | %6s %6s\n', '', 'grid', 'Newton', 'fmincon', 'HOI m/s', ...
    'hp km', 'FPA deg', 'TOF d', 'dtau', 'dthS deg', 'evals', 's');
for n = 1:numel(seeds)
    seed = seeds{n};
    [flyby, rmin_km, tmin_d] = moon_flyby(seed, halo, prm);
    if flyby
        fprintf('4. T%d screened out: lunar flyby %.0f km from the Moon, %.1f days before arrival\n', n, rmin_km, tmin_d);
        continue
    end
    tic;
    P  = dc_setup(seed.p, halo, prm);
    z0 = [0; 0; 0; -seed.tof_d*86400/prm.tstar; seed.tau; seed.thS_deg*pi/180];
    opt_newton = dc_newton(P, z0, P.newton_vars);
    if ~opt_newton.converged
        fprintf('4. T%d Newton did not converge (|F| = %.1e): seed perigee %.1f km at %.1f days\n', n, norm(opt_newton.F), seed.hmin_km, seed.tof_d);
        continue
    end
    opt_fmincon = dc_optimize(P, opt_newton);
    sec = toc;
    if opt_fmincon.converged, kept = opt_fmincon; else, kept = opt_newton; end
    z  = kept.z;
    xf = arrival_states(halo, z(5), seed.p, prm) - [0; 0; 0; z(1:3)];
    xp = deval(ode113(@(t, X) eom_bcr4bp6(t, X, prm, z(6)), [0 z(4)], xf, o6), z(4));  % independent: 6-state, tighter tolerance
    d  = xp(1:3) - [-prm.mu; 0; 0];
    fprintf('4. T%d  %9.5f %9.5f %9.5f %11.4f | %10.4f %8.1e %10.3f | %10.5f %9.4f | %6d %6.1f\n', n, seed.dv_kms, ...
        opt_newton.info.dv_kms, opt_fmincon.info.dv_kms, kept.info.hoi_kms*1e3, (norm(d) - prm.RE)*prm.lstar, asind(d.'*xp(4:6)/(norm(d)*norm(xp(4:6)))), ...
        kept.info.tof_d, z(5) - seed.tau, (z(6) - seed.thS_deg*pi/180)*180/pi, opt_fmincon.evals, sec);
end
fprintf('   grid = energy-patch TLI of the seed; Newton and fmincon totals are TLI + HOI (km/s); hp and FPA from the re-propagation\n');

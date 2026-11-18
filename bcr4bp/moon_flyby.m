function [flyby, rmin_km, tmin_d] = moon_flyby(seed, halo, prm)
% Screen a grid-search seed for a lunar flyby: a close approach to the Moon inside prm.dc.r_flyby_km on its backward arc to the seed perigee
% A flyby makes the arc too sensitive for single shooting (dc_newton, dc_optimize), so such seeds are removed rather than corrected
% Inputs: seed - struct: thS_deg, tau, p, tof_d (see correction_main); halo - halo struct (see arrival_states); prm - struct from bcr4bp_params
% Outputs: flyby - true if a close approach is inside prm.dc.r_flyby_km; rmin_km - closest approach to the Moon's center (km, Inf if none); tmin_d - its time (days before arrival)
    x0   = arrival_states(halo, seed.tau, seed.p, prm);
    f    = @(t, X) eom_bcr4bp6(t, X, prm, seed.thS_deg*pi/180);
    opts = odeset('RelTol', prm.tol, 'AbsTol', prm.tol, 'Events', @(t, X) moon_apse(t, X, prm));
    [~, ~, te, Xe] = ode113(f, [0 -seed.tof_d*86400/prm.tstar], x0, opts);

    % Close approaches: distance to the Moon has a minimum (second derivative of |r - rM|^2 positive)
    rM   = [1 - prm.mu, 0, 0];
    keep = false(size(te));
    for k = 1:numel(te)
        Xd = f(te(k), Xe(k, :).');
        keep(k) = Xe(k, 4:6)*Xe(k, 4:6).' + (Xe(k, 1:3) - rM)*Xd(4:6) > 0;
    end
    [rmin_km, tmin_d] = deal(Inf, NaN);
    if any(keep)
        r = vecnorm(Xe(keep, 1:3) - rM, 2, 2)*prm.lstar;  t = te(keep);
        [rmin_km, k] = min(r);
        tmin_d = -t(k)*prm.tstar/86400;
    end
    flyby = rmin_km < prm.dc.r_flyby_km;
end

function [g, term, dirn] = moon_apse(~, X, prm)
% Event: closest or farthest point from the Moon (logged, not terminal)
    g    = (X(1:3) - [1 - prm.mu; 0; 0]).'*X(4:6);
    term = 0;
    dirn = 0;
end

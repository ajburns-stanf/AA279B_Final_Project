function out = propagate_grid_point(x0, thS0, prm, dt_d, tol)
% Propagate one arrival state backward in the BCR4BP and record its perigees
% Inputs: x0 - arrival state from arrival_states (6x1); thS0 - Sun angle at halo arrival (rad); prm - struct from bcr4bp_params; dt_d - backward duration in days (optional, default prm.dt_d); tol - ode113 tolerance (optional, default prm.tol)
% Outputs: out - struct: tp, hp_km, Xp (perigees in backward order, Earth impact counted as altitude 0); kmin, hmin_km, tof_d, C3 (at the lowest perigee), dv_kms (burn from the circular orbit at hLEO_km to that C3); term ('none', 'earth', 'moon', 'escape'); tend
    if nargin < 4, dt_d = prm.dt_d; end
    if nargin < 5, tol = prm.tol; end
    f    = @(t, X) eom_bcr4bp6(t, X, prm, thS0);
    opts = odeset('RelTol', tol, 'AbsTol', tol, 'Events', @(t, X) grid_events(t, X, prm));
    [t, ~, te, Xe, ie] = ode113(f, [0 -dt_d*86400/prm.tstar], x0, opts);

    % True perigees: distance to Earth has a minimum (second derivative of |r - rE|^2 positive)
    keep = false(size(ie));
    for k = 1:numel(ie)
        if ie(k) == 1
            Xd = f(te(k), Xe(k, :).');
            keep(k) = Xe(k, 4:6)*Xe(k, 4:6).' + (Xe(k, 1:3) - [-prm.mu 0 0])*Xd(4:6) > 0;
        else
            keep(k) = ie(k) == 2;  % Earth impact
        end
    end
    names    = {'none', 'earth', 'moon', 'escape'};
    out.term = names{1};
    if ~isempty(ie) && ie(end) > 1, out.term = names{ie(end)}; end
    out.tend = t(end);

    out.tp = te(keep).';
    out.Xp = Xe(keep, :).';
    out.hp_km = (vecnorm(out.Xp(1:3, :) - [-prm.mu; 0; 0]) - prm.RE)*prm.lstar;
    out.hp_km(ie(keep) == 2) = 0;

    [out.hmin_km, out.kmin, out.tof_d, out.C3, out.dv_kms] = deal(Inf, 0, NaN, NaN, NaN);
    if isempty(out.tp), return; end
    [out.hmin_km, out.kmin] = min(out.hp_km);
    [out.C3, out.dv_kms] = perigee_energy(out.Xp(:, out.kmin), prm);
    out.tof_d = -out.tp(out.kmin)*prm.tstar/86400;
end

function [g, term, dirn] = grid_events(~, X, prm)
% Events: perigee or apogee (logged), Earth impact, Moon impact, escape (stop)
% Inputs: X - state (6x1), prm - struct from bcr4bp_params
% Outputs: g - event values, term - stop flags, dirn - crossing directions
    d    = X(1:3) - [-prm.mu; 0; 0];
    g    = [d.'*X(4:6); norm(d) - prm.RE; norm(X(1:3) - [1 - prm.mu; 0; 0]) - prm.RM; norm(X(1:3)) - prm.r_esc];
    term = [0; 1; 1; 1];
    dirn = [0; 0; 0; 0];
end

function [C3, dv] = perigee_energy(X, prm)
% Earth-relative launch energy and the tangential burn from the circular LEO that gives it
% Inputs: X - state at perigee (6x1), prm - struct from bcr4bp_params
% Outputs: C3 - km^2/s^2, dv - km/s (from the circular orbit at altitude hLEO_km)
    d  = X(1:3) - [-prm.mu; 0; 0];
    vI = norm(X(4:6) + cross([0; 0; 1], d))*prm.vstar;  % inertial speed relative to Earth
    C3 = vI^2 - 2*prm.GM_E/(norm(d)*prm.lstar);
    rL = prm.rLEO*prm.lstar;
    dv = sqrt(C3 + 2*prm.GM_E/rL) - sqrt(prm.GM_E/rL);
end

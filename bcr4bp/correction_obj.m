function [F, dF, J, dJ, info] = correction_obj(z, p, halo, prm)
% Earth-end residuals and total-dV cost of a two-impulse transfer, with their gradients
% The arrival state at halo phase tau (branch p), less the halo-insertion burn, is propagated backward in the BCR4BP for time tp to the departure
% Inputs: z - [dv (HOI, 3x1, VU); tp (TU, negative); tau (halo phase); thS0 (Sun angle at arrival, rad)]; p - branch (+1 exterior, -1 interior); halo - halo struct (see arrival_states); prm - struct from bcr4bp_params
% Outputs: F - [(|d| - rLEO)/rLEO; dhat.vhat] at departure (2x1, zero on the LEO with a tangential velocity); dF - dF/dz (2x6);
%          J - TLI + |HOI| (VU, HOI norm smoothed by prm.dc.hoi_eps_mps); dJ - dJ/dz (1x6);
%          info - struct: xf (arrival state), xp (departure state), tli_kms, hoi_kms, dv_kms (TLI + unsmoothed |HOI|), hp_km, fpa_deg, tof_d, thS_dep (Sun angle at departure, rad)
%          All outputs are NaN if the backward arc hits the Earth before tp
    z  = z(:);  dv = z(1:3);  tp = z(4);  tau = z(5);  thS0 = z(6);
    [xf, dxf] = arrival_states(halo, tau, p, prm);

    w0   = [xf - [0; 0; 0; dv]; reshape(eye(6), [], 1); zeros(6, 1)];  % [state; STM; dX/dthS0]
    opts = odeset('RelTol', prm.dc.tol, 'AbsTol', prm.dc.tol, 'Events', @(t, w) earth_impact(t, w, prm));
    sol  = ode113(@(t, w) eom_bcr4bp6(t, w, prm, thS0), [0 tp], w0, opts);
    info = struct('xf', xf, 'xp', nan(6, 1), 'tli_kms', NaN, 'hoi_kms', norm(dv)*prm.vstar, 'dv_kms', NaN, ...
                  'hp_km', NaN, 'fpa_deg', NaN, 'tof_d', -tp*prm.tstar/86400, 'thS_dep', thS0 + prm.thdot*tp);
    if ~isempty(sol.ie)
        [F, dF, J, dJ] = deal(nan(2, 1), nan(2, 6), NaN, nan(1, 6));
        return
    end

    w   = sol.y(:, end);
    xp  = w(1:6);  
    Phi = reshape(w(7:42), 6, 6);  
    psi = w(43:48);
    
    fp  = eom_bcr4bp6(tp, xp, prm, thS0);
    G   = [-Phi(:, 4:6), fp, Phi*dxf, psi];  % d(xp)/dz

    % Residuals: departure on the LEO sphere with zero flight-path angle (the angle is the same in the rotating and inertial frames)
    d  = xp(1:3) - [-prm.mu; 0; 0];  v = xp(4:6);
    r  = norm(d);  dh = d/r;  vn = norm(v);  vh = v/vn;
    s  = dh.'*vh;
    F  = [(r - prm.rLEO)/prm.rLEO; s];
    dF = [dh.'/prm.rLEO, zeros(1, 3); ((vh - s*dh)/r).', ((dh - s*vh)/vn).']*G;

    % Cost: tangential TLI from the circular LEO plus the HOI
    vI  = v + [-d(2); d(1); 0];  vIn = norm(vI);  % inertial velocity relative to Earth, v + z x d
    tli = vIn - sqrt((1 - prm.mu)/prm.rLEO);
    eh  = prm.dc.hoi_eps_mps/1e3/prm.vstar;
    hoi = sqrt(dv.'*dv + eh^2);
    J   = tli + hoi;
    dJ  = [(vI/vIn).'*[0 -1 0; 1 0 0; 0 0 0], (vI/vIn).']*G + [dv.'/hoi, 0, 0, 0];

    info.xp      = xp;
    info.tli_kms = tli*prm.vstar;
    info.dv_kms  = info.tli_kms + info.hoi_kms;
    info.hp_km   = (r - prm.RE)*prm.lstar;
    info.fpa_deg = asind(s);
end

function [g, term, dirn] = earth_impact(~, w, prm)
% Event: the arc reaches the Earth's surface
    g    = norm(w(1:3) - [-prm.mu; 0; 0]) - prm.RE;
    term = 1;
    dirn = 0;
end

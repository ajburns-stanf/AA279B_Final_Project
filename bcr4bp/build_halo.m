function halo = build_halo(prm, X0)
% L2 halo orbit and its stable eigenvector in the CR3BP, as in Part 1 of AA279B_FinalProject.m
% Differentially corrects x0 and vy0 at fixed z0 so that xdot = zdot = 0 at the half-period y = 0 crossing, then propagates the STM over one period
% Inputs: prm - struct from bcr4bp_params; X0 - initial guess [x; 0; z; 0; vy; 0] (optional, default the Earth-side IC of the main script)
% Outputs: halo - struct: X0 (corrected IC), P (period), sol (ode113 solution of [state; STM(:)] over one period), vs (stable eigenvector at t = 0, unit position part), lam (monodromy eigenvalues)
    if nargin < 2, X0 = [1.0606924537198088; 0; 0.07354804500048147; 0; 0.3422130947784463; 0]; end
    p0   = prm;  p0.muS = 0;  % muS = 0 is the Earth-Moon CR3BP
    f    = @(t, w) eom_bcr4bp6(t, w, p0, 0);
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    ev   = odeset(opts, 'Events', @(t, w) cross_y(t, w, -sign(X0(5))));

    x0 = X0(1);  z0 = X0(3);  vy0 = X0(5);
    for k = 1:50
        [~, ~, te, we] = ode113(f, [0 10], [x0; 0; z0; 0; vy0; 0; reshape(eye(6), [], 1)], ev);
        wf  = we(1, :).';
        Phi = reshape(wf(7:42), 6, 6);
        err = [wf(4); wf(6)];
        if norm(err) < 1e-12, break; end
        fd  = f(0, wf(1:6));
        Jh  = Phi([4 6], [1 5]) - (fd([4 6])/wf(5))*Phi(2, [1 5]);
        dX  = -Jh\err;
        x0  = x0 + dX(1);  vy0 = vy0 + dX(2);
    end
    if norm(err) >= 1e-12, error('build_halo: halo correction did not converge (|err| = %.2e)', norm(err)); end

    halo.X0  = [x0; 0; z0; 0; vy0; 0];
    halo.P   = 2*te(1);
    halo.sol = ode113(f, [0 halo.P], [halo.X0; reshape(eye(6), [], 1)], opts);
    M = reshape(halo.sol.y(7:42, end), 6, 6);  % monodromy matrix
    [V, L] = eig(M);
    halo.lam = diag(L);
    [~, is]  = min(abs(halo.lam));             % stable eigenvalue
    vs       = real(V(:, is));
    halo.vs  = vs/norm(vs(1:3));
end

function [value, isterminal, direction] = cross_y(~, w, dirn)
% Event: y = 0 crossing in direction dirn (the half period)
    value = w(2);  isterminal = 1;  direction = dirn;
end

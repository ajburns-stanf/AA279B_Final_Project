function [x0, dx0] = arrival_states(halo, tau, p, prm, eps_km)
% Arrival states near the stable manifold: halo points at phase tau shifted along the transported stable eigenvector
% Inputs: halo - struct with sol (ode113 solution of [state; STM(:)] over one period), vs (stable eigenvector at t = 0), P (period); tau - phases in [0, 1) (1xN); p - branch (+1 exterior, -1 interior); prm - struct from bcr4bp_params; eps_km - shift in km (optional, default prm.eps_km)
% Outputs: x0 - arrival states (6xN), the initial conditions for backward propagation; dx0 - derivative of x0 with respect to tau (6xN, optional)
    if nargin < 5, eps_km = prm.eps_km; end
    vs = halo.vs;
    if ~is_exterior(halo, vs, prm), vs = -vs; end  % eig sign is arbitrary

    W   = deval(halo.sol, mod(tau(:).', 1)*halo.P);
    Phi = reshape(W(7:42, :), 6, 6, []);
    vk  = reshape(pagemtimes(Phi, vs), 6, []);     % eigenvector carried to each phase

    x0 = W(1:6, :) + p*eps_km/prm.lstar*vk./vecnorm(vk(1:3, :));
    if nargout < 2, return; end

    % d/dtau = P*d/dt along the CR3BP halo: d(state)/dt = f, d(vk)/dt = A*vk
    p0  = prm;  p0.muS = 0;
    dx0 = zeros(size(x0));
    for k = 1:size(x0, 2)
        [f, A] = eom_bcr4bp6(0, W(1:6, k), p0, 0);
        n   = norm(vk(1:3, k));
        dvk = A*vk(:, k);
        du  = dvk/n - vk(:, k)*(vk(1:3, k).'*dvk(1:3))/n^3;  % derivative of vk/|vk(1:3)|
        dx0(:, k) = halo.P*(f + p*eps_km/prm.lstar*du);
    end
end

function ext = is_exterior(halo, vs, prm)
% Check if shifting along +vs at tau = 0 gives the exterior branch
% Inputs: halo - halo struct, vs - stable eigenvector at t = 0, prm - struct from bcr4bp_params
% Outputs: ext - true if the arc passes the halo's far side before the Moon's x (backward in CR3BP)
    p0 = prm;  p0.muS = 0;
    x0 = halo.sol.y(1:6, 1) + prm.eps_km/prm.lstar*vs/norm(vs(1:3));
    [~, Y] = ode113(@(t, x) eom_bcr4bp6(t, x, p0, 0), [0 -3*halo.P], x0, odeset('RelTol', 1e-10, 'AbsTol', 1e-10));
    kE  = find(Y(:, 1) > max(halo.sol.y(1, :)) + 0.05, 1);  % 0.05 safety margin past the halo's far side
    kI  = find(Y(:, 1) < 1 - prm.mu, 1);
    ext = ~isempty(kE) && (isempty(kI) || kE < kI);
end

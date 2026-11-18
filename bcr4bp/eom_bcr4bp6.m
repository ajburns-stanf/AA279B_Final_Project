function [Xdot, A] = eom_bcr4bp6(t, X, prm, thS0)
% Time derivative of BCR4BP state in the Earth-Moon rotating frame, optionally with the STM and the Sun-angle sensitivity
% Inputs: t - time; X - state [x; y; z; vx; vy; vz] (6x1), or [state; STM(:)] (42x1), or [state; STM(:); dX/dthS0] (48x1); prm - struct from bcr4bp_params; thS0 - Sun angle at t = 0 (rad)
% Outputs: Xdot - time derivative of X (same size as X); A - Jacobian of the state derivative (6x6, computed when requested or when X carries the STM)
    mu = prm.mu;  muS = prm.muS;
    x  = X(1); y  = X(2); z  = X(3);
    vx = X(4); vy = X(5); vz = X(6);

    wz    = 1;          % Earth-Moon frame rate
    r_em  = 1;          % Earth-Moon distance
    r_b1s = prm.aS;     % Sun-B1 distance
    thS   = thS0 + prm.thdot*t;
    xs = r_b1s*cos(thS);  ys = r_b1s*sin(thS);  zs = 0;

    r_ec = sqrt((x + mu)^2 + y^2 + z^2);
    r_mc = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    r_sc = sqrt((x - xs)^2 + (y - ys)^2 + (z - zs)^2);
    kS   = muS/(r_b1s/r_em)^3;

    % -2*cross(omega, v) + grad(Upsilon)
    ax =  2*wz*vy + wz^2*x - (1 - mu)*(x + mu)/r_ec^3 - mu*(x - 1 + mu)/r_mc^3 - muS*(x - xs)/r_sc^3 - kS*xs;
    ay = -2*wz*vx + wz^2*y - (1 - mu)*y/r_ec^3        - mu*y/r_mc^3            - muS*(y - ys)/r_sc^3 - kS*ys;
    az =                   - (1 - mu)*z/r_ec^3        - mu*z/r_mc^3            - muS*(z - zs)/r_sc^3 - kS*zs;

    Xdot = [vx; vy; vz; ax; ay; az];
    if numel(X) == 6 && nargout < 2, return; end

    I3 = eye(3);
    dE = [x + mu; y; z];  dM = [x - 1 + mu; y; z];  dS = [x - xs; y - ys; z - zs];
    G  = -(1 - mu)*(I3/r_ec^3 - 3*(dE*dE.')/r_ec^5) ...
         - mu*(I3/r_mc^3 - 3*(dM*dM.')/r_mc^5) ...
         - muS*(I3/r_sc^3 - 3*(dS*dS.')/r_sc^5) ...
         + diag([wz^2, wz^2, 0]);

    A  = [zeros(3), I3; G, [0 2*wz 0; -2*wz 0 0; 0 0 0]];

    if numel(X) >= 42
        Phi  = reshape(X(7:42), 6, 6);
        Xdot = [Xdot; reshape(A*Phi, [], 1)];
    end

    if numel(X) == 48
        % Sun-angle d(accel)/d(thS0) = d(accel)/d(rS) * d(rS)/d(thS)
        dadth = (GS - kS*I3)*(r_b1s*[-sin(thS); cos(thS); 0]);
        Xdot  = [Xdot; A*X(43:48) + [0; 0; 0; dadth]];
    end
end

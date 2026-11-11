function Xdot = eom_bcr4bp6(t, X, prm, thS0)
% Time derivative of BCR4BP state (no STM) in the Earth-Moon rotating frame
% Inputs: t - time, X - state [x; y; z; vx; vy; vz], prm - struct from bcr4bp_params, thS0 - Sun angle at t = 0 (rad)
% Outputs: Xdot - time derivative of X (6x1)
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
end

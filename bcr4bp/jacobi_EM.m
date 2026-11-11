function C = jacobi_EM(X, prm)
% Earth-Moon Jacobi constant from the CR3BP part of the pseudo-potential
% Inputs: X - states (Nx6, rows [x y z vx vy vz], or one 6x1 state), prm - struct from bcr4bp_params
% Outputs: C - Jacobi constant per state (Nx1)
    if iscolumn(X), X = X.'; end
    mu = prm.mu;
    x  = X(:,1); y = X(:,2); z = X(:,3);
    v2 = sum(X(:,4:6).^2, 2);

    wz   = 1;  % Earth-Moon frame rate
    r_ec = sqrt((x + mu).^2 + y.^2 + z.^2);
    r_mc = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);

    Upsilon_CR3BP = 0.5*wz^2*(x.^2 + y.^2) + (1 - mu)./r_ec + mu./r_mc;
    C = 2*Upsilon_CR3BP - v2;
end

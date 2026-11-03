%% AA 279B Final Project: L2 Halo to Earth Trajectories
% = Alexander Burns, Kaila Coimbra, Yuji Takubo =
clear; clc; close all;


%% ========= Part 1a ===============
%   Compute L2 Halo Orbit
% =================================
disp('Part 1a: Compute L2 Halo Orbit')

% -- Define System Constants --
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
mu = 0.012150585609624; % mu of the moon
C = 3.08; % jacobi constant

% -- Determine A Matrix --
a_matrix_fxn = compute_a(); % function from partialF/partialX
Phi0 = eye(6,6);

% -- Determine Liberation Points --
Ux  = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x-1+mu)./abs(x-1+mu).^3;
xL2 = fzero(Ux, [1-mu+1e-6, 1.5]);

% -- Differentially correct L2 Halo Orbit (fixed z0) --
max_iter = 50;
X0_L2_halo = [1.0606924537198088; 0; 0.07354804500048147; 0; 0.3422130947784463; 0]; % Earth-side IC
z0_h = X0_L2_halo(3);
[x0_h, vy0_h, P_L2_halo] = correct_halo(X0_L2_halo(1), z0_h, X0_L2_halo(5), mu, a_matrix_fxn, max_iter);

% -- Propagate State Transition Matrix Forward --
X0_L2_halo = [x0_h; 0; z0_h; 0; vy0_h; 0];
[t2, W2] = ode113(@(t,w) dynamics(t,w,a_matrix_fxn,mu), [0 P_L2_halo], [X0_L2_halo; Phi0(:)], options);

% -- L2 Halo Orbit Plot --
LU = 384400;                         % km (Pavlak)
mx = @(x) (x - (1-mu))*LU/1e4;       % Moon-centered, x 10^4 km
my = @(y) y*LU/1e4;
mz = @(z) z*LU/1e4;

figure(1); clf; hold on; grid on; axis equal
plot3(mx(W2(:,1)), my(W2(:,2)), mz(W2(:,3)), 'r', 'LineWidth', 1.2)  % L2 halo
plot3(0, 0, 0, 'o', 'Color',[.5 .5 .5], 'MarkerFaceColor',[.5 .5 .5])  % Moon
plot3(mx(xL2), 0, 0, 'k.', 'MarkerSize', 14)
text(mx(xL2), 0, -0.9, 'L_2'); text(0, 0, -0.9, 'Moon')
xlabel('x (10^4 km)'); ylabel('y (10^4 km)'); zlabel('z (10^4 km)')
view(3)


%% ========= Part 1b ===============
%   Compute L2 Stable Invariant Manifolds
% =================================
disp('Part 1a: Compute Stable Invariant Manifold')
eps_r = 1e-6; eps_v = 1e-6; % pos, vel perturbations

% -- Monodromy Matrix + Stable Eigenvector --
M_L2 = reshape(W2(end,7:42), 6, 6); % Monodromy matrix only valid at t = one period, Pavlak Eq. (3.69)
[V_M2, Lam_L2] = eig(M_L2); 

[~,is] = min(abs(diag(Lam_L2))); % stable (L2)
vs_L2  = real(V_M2(:,is));  vs_L2 = vs_L2/norm(vs_L2(1:3));

% -- Eigenvalue Plot --
plot_eig(Lam_L2, "L2 Halo Orbit Eigenvalues");

% -- Propagate Stable Manifold (backward in time) --
n_arcs = 100; % number of arcs to draw
time = -2.7*P_L2_halo; % propagate stable manifold backward in time for X periods
[arcsS, secS] = manifold_arcs(W2, vs_L2, eps_r, eps_v, time, a_matrix_fxn, mu, n_arcs); % only returns arcs which cross the moon's plane

% -- Plot The Arcs --
figure(3); clf; hold on; axis equal; grid on
plot3(mx(W2(:,1)), my(W2(:,2)), mz(W2(:,3)), 'k', 'LineWidth', 1.5)            % L2 halo
for i = 1:numel(arcsS)
    hS = plot3(mx(arcsS{i}(:,2)), my(arcsS{i}(:,3)), mz(arcsS{i}(:,4)), 'b');  % stable manifold
end

% Add Earth, Moon and L2 to plot
hE = plot3(mx(-mu), 0, 0, 'o', 'Color',[0 .45 .74], 'MarkerFaceColor',[0 .45 .74], 'MarkerSize', 10);  % Earth
hM = plot3(0, 0, 0, 'o', 'Color',[.4 .4 .4], 'MarkerFaceColor',[.5 .5 .5], 'MarkerSize', 6);           % Moon
hL2 = plot3(mx(xL2), 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 4);                              % L2
text(mx(-mu), 0, -2, 'Earth', 'HorizontalAlignment','center')
text(0,       0, -2, 'Moon',  'HorizontalAlignment','center')
text(mx(xL2), 0, -2, 'L_2',   'HorizontalAlignment','center')

xlabel('x (10^4 km)'); ylabel('y (10^4 km)'); zlabel('z (10^4 km)')
legend([hS hE hM hL2], {'Stable manifold','Earth','Moon','L2'}, 'Location','northeast')
view(3)


%% ========= Part 2 ===============
%   Grid Search To Identify Candidate Transfer Trajectories
% =================================


%% ======= Helper Functions =======

% ---- Part 1a: Compute A(t) (Phi_dot = A(t)*phi_0 ----
function a_matrix = compute_a()
    syms x y z vx vy vz mu real

    % -- Distances, Eqs. 2.23 - 2.24 --
    d13 = sqrt((x+mu)^2 + y^2 + z^2);
    d23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    % -- Accelerations, Eqs. 2.20 - 2.22 --
    ax =  2*vy + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
    ay = -2*vx + y - (1-mu)*y/d13^3 - mu*y/d23^3;
    az = -(1-mu)*z/d13^3 - mu*z/d23^3;

    % -- State vector and first-order EOM, Eq. 3.13 / 3.15 --
    state = [x; y; z; vx; vy; vz];
    f = [vx; vy; vz; ax; ay; az];

    % -- Jacobian, A(t) = df/dstate (Eq. 3.17) --
    A_sym = jacobian(f, state);

    % -- Optional: simplify and sanity-check --
    A_sym = simplify(A_sym);
    disp(A_sym)

    % -- Return as a function handle --
    a_matrix = matlabFunction(A_sym, 'Vars', {x, y, z, mu});
end


% ---- Part 1a: Compute Trajectory and STM ----
function eom_Aug = dynamics(t, w, a_matrix_fxn, mu)
    % Unpack w = [x; phi]' (42x1)
    x  = w(1);  y  = w(2);  z  = w(3);
    vx = w(4);  vy = w(5);  vz = w(6);
    Phi = reshape(w(7:42), 6, 6);

    % Compute x_dot Eq 2.20-2.22
    d13 = sqrt((x+mu)^2 + y^2 + z^2);
    d23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    ax =  2*vy + x - (1-mu)*(x+mu)/d13^3 - mu*(x-1+mu)/d23^3;
    ay = -2*vx + y - (1-mu)*y/d13^3 - mu*y/d23^3;
    az = -(1-mu)*z/d13^3 - mu*z/d23^3;

    x_dot = [vx; vy; vz; ax; ay; az];
    
    % Compute phi_dot = A(t)*Phi(t) - EQ 3.24
    A = a_matrix_fxn(x, y, z, mu);

    phi_dot = A*Phi;

    % Ending Augmented matrix (state + stm)
    eom_Aug = [x_dot; phi_dot(:)]; 
end


% ---- Part 1a: Differentially Correct L2 Halo Orbit at z0 ----
% Unknown: x0/vy0. Target: xdot/zdot = 0 at y = 0 crossing
% Uses STM entries (Eq. 3.24-3.25) for the Newton update.  Returns full period P.
function [x0, vy0, P] = correct_halo(x0, z0, vy0, mu, a_fxn, max_iter)
    % Define Tolerances & Events to record (crossing y = 0)
    options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @(t,w) cross_y(t,w,-sign(vy0)));
    
    % Differentially correct up to max_iter iterations
    for k = 1:max_iter
        % -- Evaluate current ICs w0 --
        w0 = [x0; 0; z0; 0; vy0; 0; reshape(eye(6),[],1)];
        [~,~,te,we,] = ode113(@(t,w) dynamics(t,w,a_fxn,mu), [0 10], w0, options);

        % -- Unpack ode113 output --
        wf  = we(1,:)'; % satellite state at the first crossing
        Phi = reshape(wf(7:42),6,6); % the corresponding 6x6 STM 

        % -- Estimate Error (Target) --
        error = [wf(4); wf(6)]; % xdot/zdot at the crossing -> which we want to go to zero

        % -- Break if Error < allowed --
        if norm(error) < 1e-12
            break
        end

        % -- If Error > allowed, perform differential correction --

        % - Compute State Derivative - 
        w_dot_full = dynamics(0, wf, a_fxn, mu); % full state derivative at crossing
        w_dot = w_dot_full(1:6); % state derivative only

        % - Compute Jacobian -
        J = Phi([4 6],[1 5]) - (w_dot([4 6])/wf(5)) * Phi(2,[1 5]);
        
        % - Update Guesses of x0 & vy0 -
        dX = -J\error;
        x0 = x0 + dX(1);
        vy0 = vy0 + dX(2);

    end
    P = 2*te(1); % symmetry -> full period = 2 * half (crossing occurs at half)
end


% ---- Part 1a: Event Function, Detect y = 0 Crossing (Half-Period) ----
function [value, isterminal, direction] = cross_y(~, w, dirn)
    value = w(2);  isterminal = 1;  direction = dirn;
end



% ---- Part 1b: Radial log map of Eigenvalues ----
function eigplot = plot_eig(lam, name)
    rho = @(r) (r>=1).*(1 + log10(max(r,1))) + (r<1)./(1 - log10(min(r,1)));
    map = @(z) z./abs(z) .* rho(abs(z));
    
    figure(); clf; hold on; box on; axis equal; title(name);
    th = linspace(0, 2*pi, 400);
    plot(cos(th), sin(th), 'Color',[.2 .2 .2], 'LineWidth', 2)   % unit circle
    xline(0, 'Color',[.5 .5 .5]); yline(0, 'Color',[.5 .5 .5])
    
    % Plot Eigenvalues (with associated color)
    for k = 1:numel(lam)
        if abs(abs(lam(k)) - 1) < 1e-3
            c = 'k';                         % on the unit circle
        elseif abs(lam(k)) > 1
            c = 'r';                         % unstable
        else
            c = 'g';                         % stable
        end
        zp = map(lam(k));                    % mapped position
        plot(real(zp), imag(zp), 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', c)
    end
    
    % -- Tick labels show the true |lambda| --
    vals = [1e-2 1e-1 1 1e1 1e2 1e3];
    pos  = rho(vals);
    xticks([-fliplr(pos) 0 pos])
    xticklabels([compose('-%g', fliplr(vals)), "0", compose('%g', vals)])
    set(gca, 'FontSize', 9)
    
    xlim([-4.5 4.5]); ylim([-2 2])
    xlabel('Re(\lambda)  (log-radial scale)'); ylabel('Im(\lambda)  (log-radial scale)')
    
    eigplot = 0;
end


% ---- Part 1b: Propagate manifold arcs ----
% Wo: orbit + STM history, vu: eigenvector at t=0, Tint<0 integrates backward (stable).
% eps_r/eps_v: position/velocity step sizes, each applied to a unit vector.
% Returns arcs (cell, [t x y z vx vy vz]) and sec rows [idx t x y z vx vy vz].
function [arcs, sec] = manifold_arcs(Wo, vu, eps_r, eps_v, Tint, a_fxn, mu, nPts)
    idx  = round(linspace(1, size(Wo,1), nPts+1));  idx(end) = [];
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    arcs = {};  sec = [];
    for j = 1:numel(idx)
        Phi_k = reshape(Wo(idx(j),7:42), 6, 6);
        vk    = Phi_k*vu;                               % transport eigenvector
        dr    = vk(1:3)/norm(vk(1:3));                  % unit position direction
        dv    = vk(4:6)/norm(vk(4:6));                  % unit velocity direction
        xk    = Wo(idx(j),1:6)';
        for s = [+1 -1]                                 % both branches
            w0 = [xk + s*[eps_r*dr; eps_v*dv]; reshape(eye(6),[],1)];
            [t,X] = ode113(@(t,w) dynamics(t,w,a_fxn,mu), [0 Tint], w0, opts);
            arcs{end+1}  = [t, X(:,1:6)];                % keep all arcs
            sec(end+1,:) = [idx(j), t(end), X(end,1:6)]; % final state of each arc
        end
    end
end
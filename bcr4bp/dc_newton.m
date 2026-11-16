function sol = dc_newton(P, z, free)
% Damped Newton on the Earth-end residual F(z) = 0, moving only the entries z(free); the other entries of z stay fixed
% The corrector uses it twice: Stage 1 (free = P.newton_vars, from the seed) and the inner solve of Stage 2 (free = P.inner, see dc_optimize)
%   Residual: F = [(|r| - rLEO)/rLEO; rhat.vhat] at departure (correction_obj), zero on the LEO with a tangential velocity (2 equations)
%   Step:     solves dF(:, free)*dz = -F in the scaled variables z(free)./scale(free). With more free variables than equations
%             (Stage 1: 4 unknowns) it is the minimum-norm step (lsqminnorm): the smallest scaled change that zeros the linearized F
%   Damping:  the step halves (up to 8 times) until |F| decreases; Newton stops if no halving does
%   Port of MDNav_agent core/guidance/newton.py
% Inputs: P - corrector setup from dc_setup (uses P.eval, P.scale, P.dc); z - starting point (6x1); free - indices of z to move
% Outputs: sol - struct: z (final point); F, dF, J, dJ, info (its evaluation, as from correction_obj); iter (Newton iterations);
%          converged (|F| < prm.dc.newton_tol)
    max_halvings = 8;
    s    = P.scale(free);
    sol  = P.eval(z);
    iter = 0;
    while norm(sol.F) >= P.dc.newton_tol && iter < P.dc.newton_iter
        iter = iter + 1;
        step = s.*lsqminnorm(sol.dF(:, free).*s.', -sol.F);
        lam  = 1;  accepted = false;
        for k = 1:max_halvings
            z_try = sol.z;  z_try(free) = sol.z(free) + lam*step;
            trial = P.eval(z_try);
            if all(isfinite(trial.F)) && norm(trial.F) < norm(sol.F)
                sol = trial;  accepted = true;
                break
            end
            lam = lam/2;
        end
        if ~accepted, break; end
    end
    sol.iter      = iter;
    sol.converged = norm(sol.F) < P.dc.newton_tol;
end

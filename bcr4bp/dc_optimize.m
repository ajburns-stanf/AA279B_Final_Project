function opt = dc_optimize(P, newton)
% Stage 2 of the differential correction (see correction_main): minimize the total dV = TLI + |HOI| over transfers that keep the Earth-end gap closed (F = 0)
%
% Idea: remove the constraint F = 0 by solving it inside the cost function (reduced-space method)
%   z = [dv (HOI, 3); tp; tau; thS0] is split into
%     outer variables, P.outer (default dv and tau): the only variables fmincon sees
%     inner variables, P.inner (default tp and thS0): recomputed from the outer ones so that F = 0
%   F = 0 is 2 equations, so fixing the 4 outer variables leaves 2 equations in the 2 inner unknowns. Each cost evaluation
%   solves them with Newton (solve_inner), so every point fmincon evaluates is a feasible transfer, and fmincon's problem is
%       min over z_o:  Jr(z_o) = J(z_o, z_i(z_o))    subject only to a box on z_o
%   The dependence z_i(z_o) still enters the gradient. Differentiating F(z_o, z_i(z_o)) = 0 gives
%       dz_i/dz_o = -F_i \ F_o,   so   dJr/dz_o = J_o + J_i*dz_i/dz_o        (F_o, F_i, J_o, J_i: derivatives w.r.t. z_o, z_i)
%
% One cost evaluation, reduced_cost(y):
%   y (scaled outer offsets) -> z_o -> solve_inner: predict z_i, then Newton until F = 0 -> Jr and dJr/dy back to fmincon
%   The inner root is unique only near the previous one (while the 2x2 F_i is nonsingular). Starting each solve from the last
%   feasible point and rejecting inner values outside their box keep it on the branch through the Stage 1 solution;
%   a rejected point returns NaN, and fmincon (sqp) shortens its step.
%
% Inputs: P - corrector setup from dc_setup; newton - Stage 1 solution from dc_newton (struct: z, F, dF, J, dJ, info, ...)
% Outputs: opt - struct: z, F, dF, J, dJ, info (the transfer at fmincon's answer, as from dc_newton); iter (fmincon iterations);
%          converged (result kept: inner solve converged, exitflag > 0 and total dV not above Newton's); exitflag;
%          evals (cost evaluations); failed (rejected points)
    outer = P.outer;  inner = P.inner;
    if numel(inner) ~= 2 || isempty(outer), error('dc_optimize: opt_inner must name 2 variables and opt_free at least one more'); end

    % Box about the Stage 1 solution: fmincon's bounds on the outer variables, and the window that accepts an inner solution
    lo = newton.z - P.box;  hi = newton.z + P.box;

    % fmincon works on y, the outer offsets from the Stage 1 solution in scaled units (one unit: 1 cm/s in dv, 1e-3 in tau)
    to_outer = @(y) newton.z(outer) + P.scale(outer).*y(:);
    y_lo = (lo(outer) - newton.z(outer))./P.scale(outer);
    y_hi = (hi(outer) - newton.z(outer))./P.scale(outer);

    % Last feasible transfer: the start of the next inner solve; updated by reduced_cost
    last = newton;
    n_evals = 0;  n_failed = 0;

    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'OptimalityTolerance', P.dc.opt_tol, ...
        'StepTolerance', 1e-10, 'MaxIterations', P.dc.opt_iter, 'MaxFunctionEvaluations', 5*P.dc.opt_iter, 'Display', P.dc.display);
    [y, ~, flag, output] = fmincon(@reduced_cost, zeros(numel(outer), 1), [], [], [], [], y_lo, y_hi, [], opts);

    % fmincon returns only y. last may be a later line-search trial, so solve the inner variables once more at y
    [opt, ok] = solve_inner(to_outer(y));
    kept = ok && flag > 0 && opt.info.dv_kms <= newton.info.dv_kms;
    if ~kept
        warning('dc_optimize:notkept', 'fmincon result not kept (exitflag %d, |F| = %.2e, dV %.4f vs Newton %.4f m/s)', ...
            flag, norm(opt.F), opt.info.dv_kms*1e3, newton.info.dv_kms*1e3);
    end
    opt.iter      = output.iterations;  % fmincon's iterations (replaces the final inner solve's)
    opt.converged = kept;
    opt.exitflag  = flag;
    opt.evals     = n_evals;
    opt.failed    = n_failed;


    function [Jr, dJr_dy] = reduced_cost(y)
    % Reduced cost: total dV (m/s) of the feasible transfer with outer variables y, and its gradient with respect to y
        n_evals = n_evals + 1;
        [sol, feasible] = solve_inner(to_outer(y));
        if ~feasible
            n_failed = n_failed + 1;
            Jr = NaN;  dJr_dy = nan(numel(outer), 1);
            return
        end
        F_o = sol.dF(:, outer);  F_i = sol.dF(:, inner);
        dzi_dzo = -(F_i\F_o);                                          % how the inner variables move to keep F = 0
        dJr_dzo = sol.dJ(outer) + sol.dJ(inner)*dzi_dzo;               % chain rule through z_i(z_o)
        Jr      = sol.J*P.vu2mps;
        dJr_dy  = (dJr_dzo.*P.scale(outer).').'*P.vu2mps;             % per scaled unit of y
        last    = sol;
    end

    function [sol, ok] = solve_inner(z_o)
    % Inner variables that close F = 0 at the outer values z_o: predictor, then Newton (corrector)
    % ok is false if Newton fails or the inner values leave their box (likely another branch of solutions)
        z = last.z;  z(outer) = z_o;
        z(inner) = last.z(inner) - last.dF(:, inner)\(last.dF(:, outer)*(z_o - last.z(outer)));  % first-order step along F = 0
        sol = dc_newton(P, z, inner);

        % check the upper/lower boundary of inner variables
        ok = sol.converged && all(sol.z(inner) >= lo(inner) & sol.z(inner) <= hi(inner));
    end
end

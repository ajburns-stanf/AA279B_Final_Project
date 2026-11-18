function P = dc_setup(p, halo, prm)
% Setup of the two-stage differential correction (dc_newton, dc_optimize) for one transfer branch
% Variables z = [dv (HOI, 3x1, VU); tp (backward time of flight, TU, negative); tau (halo phase); thS0 (Sun angle at arrival, rad)]
% Inputs: p - branch (+1 exterior, -1 interior); halo - halo struct (see arrival_states); prm - struct from bcr4bp_params
% Outputs: P - struct:
%   eval        - handle z -> struct z, F, dF, J, dJ, info: the point, its residual, cost and their gradients (correction_obj)
%   scale       - one scaled unit per entry of z (6x1, from prm.dc.scale), the units of the Newton steps and of fmincon's variables
%   box         - half-width of the box about the Newton solution per entry of z (6x1, from prm.dc.bound)
%   newton_vars - Stage 1 variables, indices into z (from prm.dc.newton_free; default dv and tp)
%   outer       - Stage 2 variables moved by fmincon (prm.dc.opt_free less prm.dc.opt_inner; default dv and tau)
%   inner       - Stage 2 variables solved by Newton at every fmincon point (prm.dc.opt_inner; default tp and thS0)
%   vu2mps      - VU to m/s (fmincon works on the cost in m/s)
%   dc          - prm.dc (tolerances and iteration limits)
    to_z = @(v) [v(1:3).'/1e3/prm.vstar; v(4)*86400/prm.tstar; v(5); v(6)*pi/180];  % prm.dc units (m/s, days, -, deg) to z units

    P.eval        = @(z) evaluate(z, p, halo, prm);
    P.scale       = to_z(prm.dc.scale);
    P.box         = to_z(prm.dc.bound);
    P.newton_vars = find(prm.dc.newton_free);
    P.inner       = prm.dc.opt_inner(:);
    P.outer       = setdiff(find(prm.dc.opt_free), P.inner);
    P.vu2mps      = 1e3*prm.vstar;
    P.dc          = prm.dc;
end

function ev = evaluate(z, p, halo, prm)
% The point z with its residual, cost and their gradients, as one struct
    ev.z = z;
    [ev.F, ev.dF, ev.J, ev.dJ, ev.info] = correction_obj(z, p, halo, prm);
end

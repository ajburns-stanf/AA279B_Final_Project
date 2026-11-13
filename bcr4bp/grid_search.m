function res = grid_search(halo, prm, th_deg, taus, dt_d, tol)
% Backward-propagate arrival states over a grid of Sun angle, halo phase and branch
% Inputs: halo - struct (see arrival_states); prm - struct from bcr4bp_params; th_deg - Sun angles at arrival in deg (1xNth); taus - halo phases (1xNtau); dt_d - backward days (optional, default prm.dt_d); tol - ode113 tolerance (optional, default prm.tol)
% Outputs: res - struct of Nth x Ntau x 2 arrays (page 1 exterior, page 2 interior): hmin_km, kmin, tof_d, C3, dv_kms, term (0 none, 1 earth, 2 moon, 3 escape), viable (hmin_km <= prm.hLEO_km); perigees (cell of [t; h_km]); settings
    if nargin < 5, dt_d = prm.dt_d; end
    if nargin < 6, tol = prm.tol; end
    nth = numel(th_deg);  ntau = numel(taus);
    X0  = cat(3, arrival_states(halo, taus, 1, prm), arrival_states(halo, taus, -1, prm));
    th  = th_deg*pi/180;

    [I, J, B] = ndgrid(1:nth, 1:ntau, 1:2);
    N = numel(I);
    [hmin, tof, C3, dv] = deal(nan(N, 1));
    [kmin, term] = deal(zeros(N, 1));
    per   = cell(N, 1);
    names = {'none', 'earth', 'moon', 'escape'};
    tic;
    parfor n = 1:N
        o = propagate_grid_point(X0(:, J(n), B(n)), th(I(n)), prm, dt_d, tol); %#ok<PFBNS> small arrays
        hmin(n) = o.hmin_km;  kmin(n) = o.kmin;  tof(n) = o.tof_d;  C3(n) = o.C3;  dv(n) = o.dv_kms;
        term(n) = find(strcmp(names, o.term)) - 1;
        per{n}  = [o.tp; o.hp_km];
    end
    sz = [nth, ntau, 2];
    res = struct('th_deg', th_deg, 'taus', taus, 'dt_d', dt_d, 'tol', tol, 'runtime_s', toc, ...
        'hmin_km', reshape(hmin, sz), 'kmin', reshape(kmin, sz), 'tof_d', reshape(tof, sz), ...
        'C3', reshape(C3, sz), 'dv_kms', reshape(dv, sz), 'term', reshape(term, sz));
    res.perigees = reshape(per, sz);
    res.viable   = res.hmin_km <= prm.hLEO_km;  % reaches LEO altitude or hits Earth (impact stored as 0 km)
end

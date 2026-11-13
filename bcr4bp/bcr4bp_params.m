function prm = bcr4bp_params(lstar)
% Return BCR4BP constants and characteristic quantities (nondimensional unless noted)
% Inputs: lstar - characteristic length in km (optional, default from the sidereal month)
% Outputs: prm - struct of constants
    GM_E  = 398600.4418;      % km^3/s^2
    GM_M  = 4902.800066;      % km^3/s^2
    GM_S  = 1.32712440018e11; % km^3/s^2
    AU    = 149597870.7;      % km
    T_sid = 27.32166155*86400; % s, sidereal month
    
    % Characteristic length lstar is Earth-Moon distance
    % Kepler's third law: circular Earth-Moon orbit with the sidereal period, so the Moon moves 1 rad per time unit
    if nargin < 1
        lstar = ((GM_E + GM_M)*(T_sid/(2*pi))^2)^(1/3); % 384747.96 km
    end

    prm.mu    = 0.012150585609624;            % JPL API value, matches halo IC
    prm.lstar = lstar;                        % km
    prm.tstar = sqrt(lstar^3/(GM_E + GM_M));  % s
    prm.vstar = prm.lstar/prm.tstar;          % km/s

    prm.muS   = GM_S/(GM_E + GM_M);               % Sun mass / Earth-Moon mass
    prm.aS    = AU/lstar;                         % Sun-B1 distance
    prm.thdot = sqrt((prm.muS + 1)/prm.aS^3) - 1; % Sun angle rate in EM frame

    prm.GM_E    = GM_E;                       % km^3/s^2
    prm.RE      = 6378.137/lstar;
    prm.RM      = 1737.4/lstar;
    prm.hLEO_km = 185;
    prm.rLEO    = (6378.137 + prm.hLEO_km)/lstar;

    % Grid-search settings
    prm.eps_km = 100;    % arrival-state shift along the stable eigenvector (km)
    prm.dt_d   = 180;    % search window (days backward)
    prm.tol    = 1e-10;  % ode113 relative and absolute tolerance
    prm.r_esc  = 8;      % escape radius from B1 (r_esc * lstar = 3.08e6 km, chosen to be twice Earth's Hill radius)
end

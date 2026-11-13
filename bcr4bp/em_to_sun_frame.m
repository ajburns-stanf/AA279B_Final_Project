function rS = em_to_sun_frame(r, thS)
% Rotate positions from the Earth-Moon rotating frame to the Sun-B1 rotating frame (origin B1, Sun on -x)
% Inputs: r - positions (Nx3, or one 3x1 position), thS - Sun angle per row (scalar or Nx1, rad)
% Outputs: rS - positions in the Sun-B1 frame (Nx3), same units as r
    if iscolumn(r), r = r.'; end
    a = pi - thS(:);  % Earth-Moon angle, Scheuerle Eq. (2.62)
    c = cos(a);  s = sin(a);
    rS = [c.*r(:,1) - s.*r(:,2), s.*r(:,1) + c.*r(:,2), r(:,3)];
end

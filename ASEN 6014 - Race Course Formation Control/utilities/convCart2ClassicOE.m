function oe = convCart2ClassicOE(cart)
% Function that converts cartesian orbital elements to their corresponding
% classical counterparts
%   Inputs:
%       - cart: Cartesian orbital element structure organized as follows:
%           - mu: Celestial body GM constant
%           - rVec: Inertial position vector
%           - vVec: Inertial velocity vector
%   Outputs:
%       - oe: Classical orbital element structure organized as follows:
%           - a: Orbit semi-major axis
%           - e: Orbit eccentricity
%           - i: Orbit inclination in radians
%           - RAAN: Orbit right ascension of the ascending node in radians
%           - argPeri: Orbit argument of periapsis in radians
%           - f: Orbit true anomaly in radians
%
%   By: Ian Faber, 08/28/2025
%

% Extract cartesian OE's
mu = cart.mu;
rVec = cart.rVec;
vVec = cart.vVec;

% Calculate radius, speed, and angular momentum
r = norm(rVec);
v = norm(vVec);

hVec = cross(rVec, vVec);
h = norm(hVec);

% Calculate semi-major axis
aInv = 2/r - (v^2)/mu;

if aInv == 0 % orbit is rectilinear/hyperbolic, a is "infinity"
    oe.a = 9e99;
else
    oe.a = 1/aInv;
end

% Calculate eccentricity
eVec = cross(vVec, hVec)/mu - rVec/r;
oe.e = norm(eVec);

% Calculate perifocal frame unit vectors and construct PN
if abs(oe.e) <= 1e-12 % We have a circular orbit
    oe.argPeri = 0;

    nVec = cross([0;0;1], hVec);
    n = norm(hVec);

    if n >= 1e-12 % not an equatorial orbit
        oe.RAAN = atan2(nVec(2), nVec(1));
    else
        oe.RAAN = 0;
    end

    oe.i = acos(hVec(3)/h);

    oe.f = atan2(rVec(2), rVec(1)); % Becomes true longitude
else
    i_e = eVec/oe.e;

    i_h = hVec/h;
    i_p = cross(i_h, i_e);
    
    PN = [i_e, i_p, i_h]';
    
    % Extract angles from PN
    oe.RAAN = atan2(PN(3,1),-PN(3,2));
    oe.i = acos(PN(3,3));
    oe.argPeri = atan2(PN(1,3), PN(2,3));
    
    % Calculate true anomaly
    i_r = rVec/r;
    oe.f = atan2(dot(cross(i_e, i_r), i_h), dot(i_e, i_r));
end


end
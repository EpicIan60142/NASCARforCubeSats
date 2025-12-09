function cart = convClassicOE2Cart(oe)
% Function that converts classical orbit elements into their cartesian
% counterparts
%   Inputs:
%       - oe: Classical orbit element structure, organized as follows:
%           - mu: Celestial body GM constant
%           - a: Orbit semi-major axis
%           - e: Orbit eccentricity
%           - i: Orbit inclination in radians
%           - RAAN: Orbit right ascension of the ascending node in radians
%           - argPeri: Orbit argument of periapsis in radians
%           - f: Orbit true anomaly in radians
%   Outputs:
%       - cart: Cartesian orbit element structure, organized as follows:
%           - rVec: Inertial position vector at provided true anomaly
%           - vVec: Inertial velocity vector at provided true anomaly
%
%   By: Ian Faber, 08/28/2025
%

% Pull out classical OE's
mu = oe.mu;
a = oe.a;
e = oe.e;
i = oe.i;
RAAN = oe.RAAN;
argPeri = oe.argPeri;
f = oe.f;

% Construct theta, r, x, y, and PN
theta = argPeri + f;

p = a*(1-e^2);
r = p/(1 + e*cos(f));

x = r*cos(f);
y = r*sin(f);

PN = EA2DCM([oe.RAAN, oe.i, oe.argPeri], [3,1,3]);

% Construct rVec
cart.rVec = PN'*[x; y; 0];

% Construct vVec
h = sqrt(mu*p);

cart.vVec = -(mu/h)*[
                        cos(RAAN)*(sin(theta) + e*sin(argPeri)) + sin(RAAN)*(cos(theta) + e*cos(argPeri))*cos(i);
                        sin(RAAN)*(sin(theta) + e*sin(argPeri)) - cos(RAAN)*(cos(theta) + e*cos(argPeri))*cos(i);
                        -(cos(theta) + e*cos(argPeri))*sin(i)
                    ];

end
function rho_H = calcLinearGenEllSol(oe_c, doe, f)
% Function that calculates the analytical solution to the linear
% generalized elliptic relative motion equations
%   Inputs:
%       - oe_c: Struct of chief spacecraft classical orbit elements as
%               calculated by convCart2ClassicOE.m. Contains the following
%               fields:
%               - a: Semi-major axis
%               - e: Eccentricity
%               - i: Inclination in radians
%               - RAAN: Right Ascension of the Ascending Node in radians
%               - argPeri: Argument of periapsis in radians
%               - f: Initial true anomaly in radians
%       - doe: Struct of orbit element differences for a deputy spacecraft
%              as follows:
%               - da: Difference in semi-major axis
%               - de: Difference in eccentricity
%               - di: Difference in inclination in radians
%               - dRAAN: Difference in RAAN in radians
%               - dargPeri: Difference in argument of periapsis in radians
%               - dM0: Difference in initial mean anomaly in radians
%       - f: True anomaly to calculate solution for in radians
%   Outputs:
%       - rho_H: Relative position vector in the Hill Frame
%
%   By: Ian Faber, 09/26/2025

    % Extract chief orbital elements
a = oe_c.a;
e = oe_c.e;
i = oe_c.i;
RAAN = oe_c.RAAN;
argPeri = oe_c.argPeri;
f0 = oe_c.f;

    % Extract orbit element differences
da = doe.da;
de = doe.de;
di = doe.di;
dRAAN = doe.dRAAN;
dargPeri = doe.dargPeri;
dM0 = doe.dM0;

    % Calculate intermediate variables for provided true anomaly and
    % elements
        % Eta
eta = sqrt(1-e^2);

        % Theta
theta = argPeri + f;

        % Eccentric and mean anomaly
E0 = convf2E(f0, e);
M0 = E0 - e*sin(E0);
E = convf2E(f, e);
M = E - e*sin(E);

dM = dM0 - (3/2)*(M - M0)*(da/a);

        % Analytical angle defs
f_u = atan2(e*dM, -eta*de);
f_v = f_u - (pi/2);
theta_w = atan2(di,-sin(i)*dRAAN);

        % Analytical delta defs
d_u = sqrt((((e^2)*(dM.^2))/(eta^2)) + de^2);
d_w = sqrt(di^2 + (sin(i)^2)*(dRAAN^2));

        % Chief orbit radius
r = (a*(1-e^2))./(1+e*cos(f));

    % Calculate nondimensional solution
u = da/a - (e*de)/(2*(eta^2)) + (d_u/(eta^2)).*(cos(f-f_u) + (e/2)*cos(2*f - f_u));
v = ((1 + (e^2)/2)*(dM/(eta^3)) + dargPeri + cos(i)*dRAAN) - (d_u/(eta^2)).*(2*sin(f - f_u) + (e/2)*sin(2*f - f_u));
w = d_w*cos(theta - theta_w);

    % Dimensionalize
x = u.*r;
y = v.*r;
z = w.*r;

    % Assign outputs
rho_H = [x; y; z];

end
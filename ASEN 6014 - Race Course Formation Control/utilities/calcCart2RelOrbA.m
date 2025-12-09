function A = calcCart2RelOrbA(oe_c, mu)
% Calculates the linear mapping matrix between cartesian and relative orbit
% element coordinates
%   Inputs:
%       - oe_c: Struct of classical orbit parameters for the chief
%               spacecraft, as calculated by convCart2ClassicOE.m. Contains
%               the following fields:
%               - a: Semi-major axis
%               - e: Eccentricity
%               - i: Inclination in radians
%               - RAAN: Right Ascension of the Ascending Node in radians
%               - argPeri: Argument of periapsis in radians
%               - f: True anomaly in radians
%       - mu: Gravitational parameter for the body of interest
%   Outputs:
%       - A: Linear mapping matrix from relative orbit elements to
%            cartesian relative coordinates
%   
%   By: Ian Faber, 09/26/2025
%

    % Extract orbit elements
a = oe_c.a;
e = oe_c.e;
i = oe_c.i;
RAAN = oe_c.RAAN;
argPeri = oe_c.argPeri;
f = oe_c.f;

    % Calculate intermediate elements
q1 = e*cos(argPeri);
q2 = e*sin(argPeri);
theta = argPeri + f;

p = a*(1-q1^2-q2^2);
r = p/(1 + q1*cos(theta) + q2*sin(theta));
h = sqrt(mu*p);

Vr = (h/p)*(q1*sin(theta) - q2*cos(theta));
Vt = (h/p)*(1 + q1*cos(theta) + q2*sin(theta));

    % Assign elements
A = zeros(6,6);

% Order:  da,            dtheta,        di,                            dq1,                                dq2,                                dRAAN
A(1,:) = [r/a,           (Vr/Vt)*r,     0,                             -(r/p)*(2*a*q1 + r*cos(theta)),     -(r/p)*(2*a*q2 + r*sin(theta)),     0                                     ]; % x row
A(2,:) = [0,             r,             0,                             0,                                  0,                                  r*cos(i)                              ]; % y row
A(3,:) = [0,             0,             r*sin(theta),                  0,                                  0,                                  -r*cos(theta)*sin(i)                  ]; % z row
A(4,:) = [-(Vr/(2*a)),   (1/r - 1/p)*h, 0,                             (1/p)*(Vr*a*q1 + h*sin(theta)),     (1/p)*(Vr*a*q2 - h*cos(theta)),     0                                     ]; % xDot row
A(5,:) = [-(3*Vt)/(2*a), -Vr,           0,                             (1/p)*(3*Vt*a*q1 + 2*h*cos(theta)), (1/p)*(3*Vt*a*q2 + 2*h*sin(theta)), Vr*cos(i)                             ]; % yDot row
A(6,:) = [0,             0,             Vt*cos(theta) + Vr*sin(theta), 0,                                  0,                                  sin(i)*(Vt*sin(theta) - Vr*cos(theta))]; % zDot row

end
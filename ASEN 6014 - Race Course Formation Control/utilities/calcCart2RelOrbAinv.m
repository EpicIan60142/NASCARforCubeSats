function Ainv = calcCart2RelOrbAinv(oe_c, mu)
% Calculates the inverse linear mapping matrix between cartesian and 
% relative orbit element coordinates, i.e. from relative cartesian to
% relative orbit elements
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
%       - Ainv: Inverse linear mapping matrix between relative orbit 
%               elements and relative cartesian coordinates, i.e. from
%               relative cartesian to relative orbit elements
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

alpha = a/r;
nu = Vr/Vt;
rho = r/p;
kappa1 = alpha*(1/rho - 1);
kappa2 = alpha*(nu^2)*(1/rho);

    % Assign elements
Ainv = zeros(6,6);

% Order:     x,                                            y,                                                                  z,                                             xDot,                      yDot,                                        zDot
Ainv(1,:) = [2*alpha*(2 + 3*kappa1 + 2*kappa2),            -2*alpha*nu*(1 + 2*kappa1 + kappa2),                                0,                                             (1/Vt)*(2*(alpha^2)*nu*p), (2*a/Vt)*(1 + 2*kappa1 + kappa2),            0                         ]; % da row
Ainv(2,:) = [0,                                            1/r,                                                                (cot(i)/r)*(cos(theta) + nu*sin(theta)),       0,                         0,                                           -(sin(theta)*cot(i))/Vt   ]; % dtheta row
Ainv(3,:) = [0,                                            0,                                                                  (sin(theta) - nu*cos(theta))/r,                0,                         0,                                           cos(theta)/Vt             ]; % di row
Ainv(4,:) = [(1/(rho*r))*(3*cos(theta) + 2*nu*sin(theta)), -(1/r)*(((nu^2)/rho)*sin(theta) + q1*sin(2*theta) - q2*cos(2*theta)), -((q2*cot(i))/r)*(cos(theta) + nu*sin(theta)), sin(theta)/(rho*Vt),       (1/(rho*Vt))*(2*cos(theta) + nu*sin(theta)), (q2*cot(i)*sin(theta))/Vt ]; % dq1 row
Ainv(5,:) = [(1/(rho*r))*(3*sin(theta) - 2*nu*cos(theta)), (1/r)*(((nu^2)/rho)*cos(theta) + q2*sin(2*theta) + q1*cos(2*theta)),  ((q1*cot(i))/r)*(cos(theta) + nu*sin(theta)),  -cos(theta)/(rho*Vt),      (1/(rho*Vt))*(2*sin(theta) - nu*cos(theta)), -(q1*cot(i)*sin(theta))/Vt]; % dq2 row
Ainv(6,:) = [0,                                            0,                                                                  -(cos(theta) + nu*sin(theta))/(r*sin(i)),      0,                         0,                                           sin(theta)/(Vt*sin(i))    ]; % dRAAN row

end
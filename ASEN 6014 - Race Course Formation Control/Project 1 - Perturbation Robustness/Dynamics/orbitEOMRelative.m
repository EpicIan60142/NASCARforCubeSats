function dX = orbitEOMRelative(t, X, X_c, intQuant, pConst, scConst, disturb)
% Function for the equations of motion for a deputy in a formation with a
% chief satellite in the presence of perturbations
%   Inputs:
%       - t: Current integration time in sec
%       - X: Deputy state vector arranged as follows:
%            [x; y; z; xDot; yDot; zDot]
%       - X_c: Chief state vector arranged as follows:
%              [X; Y; Z; Xdot; Ydot; Zdot]
%       - intQuant: Structure of intermediate quantities with the following
%                   fields:
%                   - r_c: Chief orbit radius at this instant of time
%                   - rDot_c: Chief orbit radius rate of change
%                   - fDot: True anomaly rate of change
%                   - HN: DCM to the chief's hill frame from the inertial
%                         frame
%       - pConst: Planetary constant vector containing mu, Ri, J2, and
%                 exponential atmosphere model parameters
%       - scConst: Spacecraft parameter structure containing mass,
%                  coefficient of drag, and cross-sectional area
%       - disturb: Boolean of disturbances to include or ignore in the
%                  following order:
%           - J2: Include (true), ignore (false)
%           - Drag: Include (true), ignore (false)
%   Outputs:
%       - dX: State rate of change vector as follows:
%             [xDot; yDot; zDot; xDDot; yDDot; zDDot]
%
%   By: Ian Faber, 08/26/2025
%

% Extract states from state vectors
x = X(1);
y = X(2);
z = X(3);
xDot = X(4);
yDot = X(5);
zDot = X(6);

X_d_N = convDeputyH2N(X_c, X, pConst);

rVec_d = X_d_N(1:3);
vVec_d = X_d_N(4:6);

cart.rVec = rVec_d;
cart.vVec = vVec_d;
cart.mu = pConst.mu;
oe_d = convCart2ClassicOE(cart);

% Extract intermediate quantities
r_c = intQuant.r_c;
rDot_c = intQuant.rDot_c;
fDot = intQuant.fDot;
HN = intQuant.HN;

% Calculate intermediate quantities
    % DCM to deputy hill frame from inertial frame
r_d = norm(rVec_d);
v_d = norm(vVec_d);

hVec_d = cross(rVec_d, vVec_d);
h_d = norm(hVec_d);

i_r_d = rVec_d/r_d;
i_h_d = hVec_d/h_d;
i_theta_d = cross(i_h_d, i_r_d);

HdN = [i_r_d'; i_theta_d'; i_h_d'];

    % DCM to chief hill frame from deputy hill frame
HHd = HN*HdN';

    % Deputy eccentricity
e_d = oe_d.e;

    % Deputy true anomaly
f_d = oe_d.f;

    % Deputy semi-latus rectum
p_d = oe_d.a*(1-e_d^2);

% Extract disturbance booleans
J2_enable = disturb(1);
if length(disturb) > 1
    drag_enable = disturb(2);
else
    drag_enable = false;
end

% Calculate unperturbed acceleration
xDDot = 2*fDot*(yDot - y*(rDot_c/r_c)) + x*fDot^2 + pConst.mu/(r_c^2) - (pConst.mu/(r_d^3))*(r_c + x);
yDDot = -2*fDot*(xDot - x*(rDot_c/r_c)) + y*fDot^2 - (pConst.mu/(r_d^3))*y;
zDDot = -(pConst.mu/(r_d^3))*z;

a_unperturb = [xDDot; yDDot; zDDot];
a_perturb = zeros(3,1);

if J2_enable
    A_J2 = [ 
                (15/2)*(rVec_d(3)/r_d)^2-(3/2)      0                              0;
                0                               (15/2)*(rVec_d(3)/r_d)^2-(3/2)     0;
                0                               0                              (15/2)*(rVec_d(3)/r_d)^2-(9/2)
           ];

    a_perturb = a_perturb + ((pConst.Ri/r_d)^2)*(pConst.mu/(r_d^3))*pConst.J2*HN*A_J2*rVec_d;
end

if drag_enable
    rho = calcDensity(pConst, r_d);

    a_perturb = a_perturb - 0.5*rho*(scConst.Cd*scConst.A/scConst.m)*v_d*HHd*[h_d*e_d*sin(f_d)/p_d; h_d/r_d; 0];
end

accel = a_unperturb + a_perturb;

dX = [xDot; yDot; zDot; accel];

end
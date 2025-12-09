function [dX, u] = ringFormationFeedbackControl(t, X, rho, kConst, pConst, scConst)
% Function that implements an inertial feedback control law for race course
% ring formation control, feeding forward for a desired ring center from
% the center of the race course
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; X_d]
%       - rho: Desired ring center position in the Hill Frame like so:
%              rho = [x; y; z]
%       - kConst: Structure of control gains with K1 and K2 as fields
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%       - scConst: Ring structure containing drag parameters Cd, m, and A
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%       - u: Control in the body frame: [u_r; u_theta; u_h]
%
%   By: Ian Faber, 11/19/2025

% Extract states
X_c = X(1:6);
X_d = X(7:12);

% Calculate intermediate variables
rChief = X_c(1:3);
vChief = X_c(4:6);
h = norm(cross(rChief, vChief));
fDot = h/(norm(rChief)^2);

rHat = rChief/norm(rChief);
hHat = cross(rChief, vChief)/h;
thetaHat = cross(hHat, rHat);

NH = [rHat, thetaHat, hHat];

r_c = norm(rChief);
rDot_c = dot(NH'*vChief, [1;0;0]); % Radius rate of change

% Calculate reference inertial state
x = rho(1);
y = rho(2);

r_r = rChief + NH*rho;
v_r = vChief -y*fDot*rHat + x*fDot*thetaHat;

X_r = [r_r; v_r];

% Calculate error vector
dX_d = X_d - X_r;

% Calculate control effort
    % Gains
K1 = kConst.K1;
K2 = kConst.K2;

    % Natural accel
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);
fDeputy = orbitEOM(t, X_d, pConst, scConst, true(1,2));
fDeputy = fDeputy(4:6);
fDesired = orbitEOM(t, X_r, pConst, scConst, true(1,2));
fDesired = fDesired(4:6);

    % Feedforward
A = [
        fDot,           2*rDot_c/r_c,   0;
        2*rDot_c/r_c,   fDot,           0;
        0,              0,              0
    ];
u_d = fChief - fDesired - fDot*NH*A*rho;

    % Control law
u = -(fDeputy - fDesired) - K1*dX_d(1:3) - K2*dX_d(4:6) + u_d;

% Assign output
dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

dX = [dX_c; dX_d];

% Report control in body frame
u = NH'*u;

end
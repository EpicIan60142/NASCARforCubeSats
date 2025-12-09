function dX = inertialFeedbackControlRadFF(t, X, offset, kConst, pConst)
% Function that implements an inertial feedback control law for relative
% motion between a chief and deputy spacecraft, feeding forward for a
% constant radial offset from the chief
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; X_d]
%       - offset: Radial offset in the radial direction
%       - kConst: Structure of control gains with K1 and K2 as fields
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%
%   By: Ian Faber, 11/06/2025

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

rDot_c = dot(NH'*vChief, [1;0;0]); % Radius rate of change

% Calculate reference inertial state
r_r = rChief + offset*rHat;
v_r = vChief + offset*fDot*thetaHat;

X_r = [r_r; v_r];

% Calculate error vector
dX_d = X_d - X_r;

% Calculate control effort
    % Gains
K1 = kConst.K1;
K2 = kConst.K2;

    % Natural accel
fDeputy = orbitEOM(t, X_d, pConst, struct(), false(1,2));
fDeputy = fDeputy(4:6);
fDesired = orbitEOM(t, X_r, pConst, struct(), false(1,2));
fDesired = fDesired(4:6);

    % Feedforward
u_d = pConst.mu*((r_r/norm(r_r)^3) - (rChief/norm(rChief)^3)) - NH*[offset*fDot^2; 2*offset*fDot*norm(rDot_c)/norm(rChief); 0];

    % Control law
u = -(fDeputy - fDesired) - K1*dX_d(1:3) - K2*dX_d(4:6) + u_d;

% Assign output
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);

dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

dX = [dX_c; dX_d];

end
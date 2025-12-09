function [dX, u, X_r_Hill] = LVLHFeedbackControlCoeffs(t, X, coeffs, kConst, pConst)
% Function that implements an LVLH feedback control law for relative
% motion between a chief and deputy spacecraft, with desired motion 
% specified with CWH analytical solution coefficients.
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; x_d_H], where x_d_H is in the hill frame
%       - coeffs: Structure of analytical CWH coefficients with the
%                 following fields: A0, B0, xOff, yOff, alpha, beta
%       - kConst: Structure of control gains with K and P as fields
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%
%   By: Ian Faber, 11/06/2025

% Extract states
X_c = X(1:6);
X_d_Hill = X(7:12);

% Calculate orbit elements of use
    % Extract chief params
mu = pConst.mu;
rVec = X_c(1:3);
vVec = X_c(4:6);
rChief = norm(rVec);

    % Perform nonlinear mapping
cart_c.mu = mu;
cart_c.rVec = rVec;
cart_c.vVec = vVec;
oe_c = convCart2ClassicOE(cart_c);

    % Calculate q1, q2, and theta
q1 = oe_c.e*cos(oe_c.argPeri);
q2 = oe_c.e*sin(oe_c.argPeri);
theta = oe_c.argPeri + oe_c.f;

    % Calculate thetaDot and thetaDDot
h = norm(cross(rVec, vVec));
thetaDot = h/(rChief^2); % Assumes argPeri is constant
thetaDDot = -2*(mu/(rChief^3))*(q1*sin(theta) - q2*cos(theta));

% Calculate dynamics matrices
    % A1
A1 = [
        2*(mu/(rChief^3)) + thetaDot^2, thetaDDot,                    0;
        -thetaDDot,                     thetaDot^2 - (mu/(rChief^3)), 0;
        0,                              0,                            -mu/(rChief^3)
     ];

    % A2
A2 = [
        0,              2*thetaDot, 0;
        -2*thetaDot,    0,          0;
        0,              0,          0
     ];

% Calculate reference inertial state
A0 = coeffs.A0;
B0 = coeffs.B0;
xOff = coeffs.xOff;
yOff = coeffs.yOff;
alpha = coeffs.alpha;
beta = coeffs.beta;

n = sqrt(pConst.mu/(oe_c.a^3));

x_r = A0*cos(n*t + alpha) + xOff;
y_r = -2*A0*sin(n*t + alpha) - 1.5*n*t*xOff + yOff;
z_r = B0*cos(n*t + beta);

xDot_r = -A0*n*sin(n*t + alpha);
yDot_r = -2*A0*n*cos(n*t + alpha) - 1.5*n*xOff;
zDot_r = -B0*n*sin(n*t + beta);

X_r_Hill = [x_r; y_r; z_r; xDot_r; yDot_r; zDot_r];

% Calculate error vector
deltaX = X_d_Hill - X_r_Hill;

% Calculate control effort
K = kConst.K;
P = kConst.P;

u = -[A1 + K*eye(size(A1)), A2 + P*eye(size(A2))]*deltaX;

% Assign output
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);

fDeputy = A1*X_d_Hill(1:3) + A2*X_d_Hill(4:6);
    % CWH equations
% fDeputy = [
%             2*n*X_d_Hill(5) + 3*(n^2)*X_d_Hill(1);
%             -2*n*X_d_Hill(4);
%             -n^2*X_d_Hill(3);
%           ];

dX_c = [X_c(4:6); fChief];
dX_d = [X_d_Hill(4:6); fDeputy + u];

dX = [dX_c; dX_d];

end
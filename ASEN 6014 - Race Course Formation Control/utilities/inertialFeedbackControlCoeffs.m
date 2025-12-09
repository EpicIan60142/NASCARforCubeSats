function dX = inertialFeedbackControlCoeffs(t, X, coeffs, kConst, pConst, distConst)
% Function that implements an inertial feedback control law for relative
% motion between a chief and deputy spacecraft, with desired motion 
% specified with CWH analytical solution coefficients.
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; X_d]
%       - coeffs: Structure of analytical CWH coefficients with the
%                 following fields: A0, B0, xOff, yOff, alpha, beta
%       - kConst: Structure of control gains with K1 and K2 as fields
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%       - distConst: Structure of disturbance coefficients with the
%                    following fields: aDragCoeff
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%
%   By: Ian Faber, 11/06/2025

% Extract states
X_c = X(1:6);
X_d = X(7:12);

% Calculate orbit elements of use
mu = pConst.mu;
rVec = X_c(1:3);
vVec = X_c(4:6);
r = norm(rVec);
v = norm(vVec);

    % Calculate semi-major axis
aInv = 2/r - (v^2)/mu;

if aInv == 0 % orbit is rectilinear/hyperbolic, a is "infinity"
    a = 9e99;
else
    a = 1/aInv;
end

% Calculate reference inertial state
A0 = coeffs.A0;
B0 = coeffs.B0;
xOff = coeffs.xOff;
yOff = coeffs.yOff;
alpha = coeffs.alpha;
beta = coeffs.beta;

n = sqrt(pConst.mu/(a^3));

x_r = A0*cos(n*t + alpha) + xOff;
y_r = -2*A0*sin(n*t + alpha) - 1.5*n*t*xOff + yOff;
z_r = B0*cos(n*t + beta);

xDot_r = -A0*n*sin(n*t + alpha);
yDot_r = -2*A0*n*cos(n*t + alpha) - 1.5*n*xOff;
zDot_r = -B0*n*sin(n*t + beta);

X_r_Hill = [x_r; y_r; z_r; xDot_r; yDot_r; zDot_r];

X_r = convDeputyH2N(X_c, X_r_Hill, pConst);

% Calculate error vector
dX_d = X_d - X_r;

% Calculate control effort
K1 = kConst.K1;
K2 = kConst.K2;

fDeputy = orbitEOM(t, X_d, pConst, struct(), false(1,2));
fDeputy = fDeputy(4:6);
fDesired = orbitEOM(t, X_r, pConst, struct(), false(1,2));
fDesired = fDesired(4:6);

u = -(fDeputy - fDesired) - K1*dX_d(1:3) - K2*dX_d(4:6);

% Assign output
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);

dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

if abs(distConst.aDragCoeff) > 0
    dX_d(4:6) = dX_d(4:6) - distConst.aDragCoeff*(X_d(4:6)/norm(X_d(4:6)));
end

dX = [dX_c; dX_d];

end
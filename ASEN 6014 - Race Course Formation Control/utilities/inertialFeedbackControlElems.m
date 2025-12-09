function dX = inertialFeedbackControlElems(t, X, doe_r, kConst, pConst)
% Function that implements an inertial feedback control law for relative
% motion between a chief and deputy spacecraft, with desired motion 
% specified with orbit element differences.
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; X_d]
%       - doe_r: Structure of desired orbit element differences with the
%                following fields: da, de, di, dRAAN, dargPeri, dM0
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

% Calculate reference elements
cart.rVec = X_c(1:3);
cart.vVec = X_c(4:6);
cart.mu = pConst.mu;
oe_c = convCart2ClassicOE(cart);

oe_r.mu = pConst.mu;
oe_r.a = oe_c.a + doe_r.da;
oe_r.e = oe_c.e + doe_r.de;
oe_r.i = oe_c.i + doe_r.di;
oe_r.RAAN = oe_c.RAAN + doe_r.dRAAN;
oe_r.argPeri = oe_c.argPeri + doe_r.dargPeri;
M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
M_r = M_c + doe_r.dM0;
E_r = convM2E(M_r, oe_r.e, false);
oe_r.f = convE2f(E_r, oe_r.e);

% Calculate reference inertial deputy state
cart_r = convClassicOE2Cart(oe_r);
X_r = [cart_r.rVec; cart_r.vVec];

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

dX = [dX_c; dX_d];

end
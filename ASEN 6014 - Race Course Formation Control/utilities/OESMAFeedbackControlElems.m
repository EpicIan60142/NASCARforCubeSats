function dX = OESMAFeedbackControlElems(t, X, doe_r, kConst, pConst)
% Function that implements a semi-major axis only feedback control law for
% relative motion between a chief and deputy spacecraft, with desired 
% motion specified with orbit element differences.
%   Inputs:
%       - t: Current integration time in sec
%       - X: Current combined state vector as follows:
%            X = [X_c; X_d]
%       - doe_r: Structure of desired orbit element differences with any of 
%                the following fields, not necessarily all: 
%                   da, de, di, dRAAN, dargPeri, dM0
%                If a field is missing, it is ignored and assumed not
%                desirable to control.
%       - kConst: Structure of control gains with the following fields:
%                 - K: semi major axis gain
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

% Calculate state elements
cart.rVec = X_c(1:3);
cart.vVec = X_c(4:6);
cart.mu = pConst.mu;
oe_c = convCart2ClassicOE(cart);

cart.rVec = X_d(1:3);
cart.vVec = X_d(4:6);
cart.mu = pConst.mu;
oe_d = convCart2ClassicOE(cart);

rDeputy = cart.rVec;
vDeputy = cart.vVec;

% Calculate DCM
rHat = rDeputy/norm(rDeputy);
hHat = cross(rDeputy, vDeputy)/norm(cross(rDeputy, vDeputy));
thetaHat = cross(hHat, rHat);

NHd = [rHat, thetaHat, hHat];

% Calculate reference elements
oe_r.a = oe_c.a + doe_r.da;

% Assign B matrix
r = norm(rDeputy);
h = norm(cross(rDeputy, vDeputy));
p = oe_d.a*(1-oe_d.e^2);
eta = sqrt(1-oe_d.e^2);

B = [
        (2*oe_d.a^2*oe_d.e*sin(oe_d.f))/(h),              (2*oe_d.a^2*p)/(h*r),                0;
        (p*sin(oe_d.f))/h,                              ((p+r)*cos(oe_d.f) + r*oe_d.e)/h,    0;
        0,                                              0,                                   (r*cos(oe_d.argPeri + oe_d.f))/h;
        0                                               0,                                   (r*sin(oe_d.argPeri + oe_d.f))/(h*sin(oe_d.i));
        -(p*cos(oe_d.f))/h,                             ((p+r)*sin(oe_d.f))/(h*oe_d.e),      -(r*sin(oe_d.argPeri + oe_d.f)*cos(oe_d.i))/(h*sin(oe_d.i));
        (eta*(p*cos(oe_d.f)-2*r*oe_d.e))/(h*oe_d.e),    -(eta*(p+r)*sin(oe_d.f))/(h*oe_d.e), 0
    ];

% Assign gains
K = kConst.K;

% Calculate control effort and convert to inertial frame
da = oe_d.a - oe_r.a;
u = -K*B(1,:)'*da;

u = NHd*u;

% Assign output
fChief = orbitEOM(t, X_c, pConst, struct(), false(1,2));
fChief = fChief(4:6);

fDeputy = orbitEOM(t, X_d, pConst, struct(), false(1,2));
fDeputy = fDeputy(4:6);

dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

dX = [dX_c; dX_d];

end
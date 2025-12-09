function [dX, u, doe, oe_d, oe_c] = OEAltFeedbackControlElems(t, X, doe_r, kConst, pConst, J2Dist)
% Function that implements an alternate orbit element feedback control law 
% for relative motion between a chief and deputy spacecraft, with desired 
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
%                 - P11, P22, P33: control gains
%       - pConst: Planetary constants vector containing mu for the
%                 celestial body of interest
%       - J2Dist: Boolean for whether to apply J2 dynamics or not
%   Outputs:
%       - dX: Rate of change vector as follows:
%             [dX_c; dX_d]
%       - u: Control effort for the given time and state
%       - doe: Controller orbit element difference vector at the given time
%       - oe_d: Deputy orbit elements for the given time and state
%       - oe_c: Chief orbit elements for the given time and state
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

rChief = cart.rVec;
vChief = cart.vVec;

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

% Calculate reference elements and populate B matrix and difference vector
B = [];
doe = [];
h = norm(cross(X_d(1:3), X_d(4:6)));
r_e = 6378; % km, radius of Earth
r = norm(X_d(1:3));
p = oe_d.a*(1-oe_d.e^2);
eta = sqrt(1-oe_d.e^2);

oe_r.mu = pConst.mu;
try doe_r.da;
    oe_r.a = oe_c.a + doe_r.da;
    Belem = [(2*oe_d.a^2*oe_d.e*sin(oe_d.f))/(h*r_e), (2*oe_d.a^2*p)/(h*r*r_e), 0];
    B = [B; Belem];
    doe = [doe; oe_d.a/r_e - oe_r.a/r_e];
end

try doe_r.de;
    oe_r.e = oe_c.e + doe_r.de;
    Belem = [(p*sin(oe_d.f))/h, ((p+r)*cos(oe_d.f) + r*oe_d.e)/h, 0];
    B = [B; Belem];
    doe = [doe; oe_d.e - oe_r.e];
end

try doe_r.di;
    oe_r.i = oe_c.i + doe_r.di;
    Belem = [0, 0, (r*cos(oe_d.argPeri + oe_d.f))/h];
    B = [B; Belem];
    doe = [doe; oe_d.i - oe_r.i];
end

try doe_r.dRAAN;
    oe_r.RAAN = oe_c.RAAN + doe_r.dRAAN;
    Belem = [0, 0, (r*sin(oe_d.argPeri + oe_d.f))/(h*sin(oe_d.i))];
    B = [B; Belem];
    doe = [doe; oe_d.RAAN - oe_r.RAAN];
end

try doe_r.dargPeri;
    oe_r.argPeri = oe_c.argPeri + doe_r.dargPeri;
    Belem = [-(p*cos(oe_d.f))/(h*oe_d.e), ((p+r)*sin(oe_d.f))/(h*oe_d.e), -(r*sin(oe_d.argPeri + oe_d.f)*cos(oe_d.i))/(h*sin(oe_d.i))];
    B = [B; Belem];
    doe = [doe; oe_d.argPeri - oe_r.argPeri];
end

try doe_r.dM0;
    M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
    M_r = M_c + doe_r.dM0;
    M_d = convE2M(convf2E(oe_d.f, oe_d.e), oe_d.e);
    Belem = [(eta*(p*cos(oe_d.f)-2*r*oe_d.e))/(h*oe_d.e), -(eta*(p+r)*sin(oe_d.f))/(h*oe_d.e), 0];
    B = [B; Belem];
    doe = [doe; M_d - M_r];
end

% Assign gains
P = diag([kConst.P11, kConst.P22, kConst.P33]);

% Calculate control effort and convert to inertial frame
if oe_d.e > 0.999
    u = zeros(3,1);
    t
else
    u = -P*B'*doe;

    u = NHd*u;
end

% Assign output
disturb = [J2Dist, false];
fChief = orbitEOM(t, X_c, pConst, struct(), disturb);
fChief = fChief(4:6);

fDeputy = orbitEOM(t, X_d, pConst, struct(), disturb);
fDeputy = fDeputy(4:6);

dX_c = [X_c(4:6); fChief];
dX_d = [X_d(4:6); fDeputy + u];

dX = [dX_c; dX_d];

end
function dX = CHWEOM(t, X, cubeSatParams, courseParams)
% Function that implements the Clohessy-Hill-Wiltshire equations, including
% the optimal control law for the CubeSat racing problem
%   Inputs:
%       - t: Current integration time
%       - X: Current CubeSat state, including adjoints:
%            [x; y; z; vx; vy; vz; px; py; pz; pvx; pvy; pvz]
%       - cubeSatParams: Parameters structure for the CubeSat of interest,
%                        as defined by generateCubesat.m
%       - courseParams: Parameters structure for the current race course
%   Outputs:
%       - dX: Rate of change of CubeSat state including adjoints
%
%   By: Ian Faber, 03/23/2025
%
%% Parse state
x = X(1); 
y = X(2); 
z = X(3); 
vx = X(4); 
vy = X(5); 
vz = X(6); 
px = X(7);
py = X(8);
pz = X(9);
pvx = X(10);
pvy = X(11);
pvz = X(12);

%% Get constants
n = courseParams.n;

uMax = 0;
if length(cubeSatParams.uMax) > 1 % Axial thrusting is at play!
    u = cubeSatParams.uMax;
    if abs(pvx) > 1 % Cubesat has max x acceleration at time t
        uMax = uMax + [u(1); 0; 0];
    end
    if abs(pvy) > 1 % Cubesat has max y acceleration at time t
        uMax = uMax + [0; u(2); 0];
    end
    if abs(pvz) > 1 % Cubesat has max z acceleration at time t
        uMax = uMax + [0; 0; u(3)];
    end
else
    uMax = cubeSatParams.uMax;
end

%% Physical dynamics
rDot = [vx; vy; vz];

pvHat = [pvx; pvy; pvz]/norm([pvx; pvy; pvz]);
g = [2*n*vy + 3*n^2*x; -2*n*vx; -n^2*z];
vDot = g - norm(uMax)*pvHat;

%% Adjoint dynamics
pr = [px; py; pz];
pv = [pvx; pvy; pvz];

prDot = -[3*n^2 0 0; 0 0 0; 0 0 -n^2]'*pv;
pvDot = -pr - [0 2*n 0; -2*n 0 0; 0 0 0]'*pv;

%% Assign output
dX = [rDot; vDot; prDot; pvDot];

end
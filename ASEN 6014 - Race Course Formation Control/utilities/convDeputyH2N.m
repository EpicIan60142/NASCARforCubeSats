function x_d_N = convDeputyH2N(x_c_N, x_d_H, const)
% Function that converts a deputy spacecraft state from Hill frame
% coordinates to inertial coordinates
%   Inputs:
%       - x_c_N: Chief spacecraft state in inertial coordinates as follows:
%                x_c_N = [r_c_N; v_c_N]
%       - x_d_H: Deputy spacecraft state in Hill frame coordinates as
%                follows:
%                x_d_H = [rho_d_H; rhoPrime_d_H]
%       - const: Constants structure containing a const.mu field
%   Outputs:
%       - x_d_N: Deputy spacecraft state in inertial frame coordinates as
%                follows:
%                x_d_N = [r_d_N; v_d_N]
%
%   By: Ian Faber, 09/14/2025
%

% % Extract orbital elements from chief spacecraft
% cart.rVec = x_c_N(1:3); 
% cart.vVec = x_c_N(4:6);
% cart.mu = const.mu;
% oe = convCart2ClassicOE(cart);
% 
% % Create DCM to Hill frame from inertial frame
% HN = EA2DCM([oe.RAAN, oe.i, oe.argPeri + oe.f], [3,1,3]);
% NH = HN';

% Extract position and velocity
r_c_N = x_c_N(1:3);
v_c_N = x_c_N(4:6);

rho_d_H = x_d_H(1:3);
rhoPrime_d_H = x_d_H(4:6);

% Compute Hill Frame
rHat = r_c_N/norm(r_c_N);
hHat = cross(r_c_N, v_c_N)/norm(cross(r_c_N, v_c_N));
thetaHat = cross(hHat, rHat);

NH = [rHat, thetaHat, hHat];

% Calculate orbital rate and angular velocity
fDot = norm(cross(r_c_N, v_c_N))/(norm(r_c_N)^2); % fDot = h/r_c^2

w_ON_H = [0; 0; fDot];

% Map vectors
r_d_N = r_c_N + NH*rho_d_H;
v_d_N = v_c_N + NH*(rhoPrime_d_H + cross(w_ON_H, rho_d_H));

% Assign output
x_d_N = [r_d_N; v_d_N];

end
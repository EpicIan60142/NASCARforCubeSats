function dX = orbitEOMInertial(t, X, pConst, scConst, disturb)
% Function for the equations of motion for a satellite in orbit around
% earth, with the option of enabling J2 perturbations
%   Inputs:
%       - t: Current integration time in sec
%       - X: State vector arranged as follows:
%            [X; Y; Z; Xdot; Ydot; Zdot]
%       - pConst: Planetary constant vector containing mu, Ri, and J2
%       - scConst: Spacecraft parameter structure containing mass,
%                  coefficient of drag, and cross-sectional area
%       - disturb: Boolean of disturbances to include or ignore in the
%                  following order:
%           - J2: Include (true), ignore (false)
%           - Drag: Include (true), ignore (false)
%   Outputs:
%       - dX: State rate of change vector as follows:
%             [Xdot; Ydot; Zdot; Xddot; Yddot; Zddot]
%
%   By: Ian Faber, 08/26/2025
%

% Extract states from state vector
rVec = X(1:3);
r = norm(rVec);

vVec = X(4:6);

% Extract disturbance booleans
J2_enable = disturb(1);
if length(disturb) > 1
    drag_enable = disturb(2);
else
    drag_enable = false;
end

% Calculate unperturbed acceleration
a_unperturb = -(pConst.mu/(r^3))*rVec;
a_perturb = zeros(3,1);

if J2_enable
    A_J2 = [ 
                (15/2)*(rVec(3)/r)^2-(3/2)      0                              0;
                0                               (15/2)*(rVec(3)/r)^2-(3/2)     0;
                0                               0                              (15/2)*(rVec(3)/r)^2-(9/2)
           ];

    a_perturb = a_perturb + ((pConst.Ri/r)^2)*pConst.J2*A_J2*(pConst.mu/(r^3))*rVec;
end

if drag_enable
    v = norm(vVec);
    rho = calcDensity(pConst, r);

    a_perturb = a_perturb - 0.5*rho*(scConst.Cd*scConst.A/scConst.m)*v*vVec;
end

accel = a_unperturb + a_perturb;

dX = [vVec; accel];

end
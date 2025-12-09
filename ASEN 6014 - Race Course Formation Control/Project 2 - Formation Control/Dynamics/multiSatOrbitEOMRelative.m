function dX = multiSatOrbitEOMRelative(t, X, pConst, scConst, disturb)
% Function that integrates chief motion and multiple deputies at the same
% time
%   Inputs:
%       - t: Current integration time
%       - X: State vector of n stacked satellite states, like so:
%            [X_chief; X_deputy1; ...; X_deputyN], where 
%            X_chief = [X; Y; Z; Xdot; Ydot; Zdot] in the inertial frame 
%              and
%            X_deputyN = [x; y; z; xDot; yDot; zDot] in the Hill frame
%       - pConst: Planetary constants structure containing Ri, mu, and J2
%       - scConst: Vector of spacecraft parameter structures containing
%                  mass, coefficient of drag, and cross-sectional area
%       - disturb: Boolean of disturbances to include or ignore in the
%                  following order for each satellite:
%           - J2: Include (true), ignore (false)
%           - Drag: Include (true), ignore (false)
%                  disturb = [disturb_sat1, disturb_sat2, ...,
%                             disturb_satn]
%   Outputs:
%       - dX: Rate of change of the state vectors, like so:
%             [dX_1; dX_2; ...; dX_n]
%
%   By: Ian Faber, 10/07/2025
%

% Number of states being integrated
nStates = 6;

% Extract chief state
X_chief = X(1:nStates);

% Extract deputies and determine the number that exist
X_deputies = X(nStates+1:end);

nDeputies = length(X_deputies)/nStates;

% Prepare rate of change vector
dX = [];

% Calculate chief rate of change and add to overall rate of change vector
dX_chief = orbitEOMInertial(t, X_chief, pConst, scConst(1), disturb(:,1));

dX = [dX; dX_chief];

% Calculate intermediate quantities
    % Calculate DCM from inertial to Hill Frame
rVec = X_chief(1:3);
vVec = X_chief(4:6);

i_r = rVec/norm(rVec);
i_h = (cross(rVec, vVec)/norm(cross(rVec,vVec)));
i_theta = cross(i_h, i_r);

HN = [i_r'; i_theta'; i_h'];

    % Calculate chief radius and its rate of change
v_c_H = HN*vVec;
r_c = norm(rVec);
rDot_c = dot(v_c_H, [1;0;0]);

    % Calculate rate of change of true anomaly
h = norm(cross(rVec, vVec));
fDot = h/(r_c^2);

    % Compile into a structure
intQuant.r_c = r_c;
intQuant.rDot_c = rDot_c;
intQuant.fDot = fDot;
intQuant.HN = HN;

% Calculate EOM for each satellite
for k = 1:nDeputies
    idx = (6*k + 1):(6*k + 6);
    X_i = X(idx);
    dX_i = orbitEOMRelative(t, X_i, X_chief, intQuant, pConst, scConst(k+1), disturb(:,k+1));
    dX = [dX; dX_i];
end

end
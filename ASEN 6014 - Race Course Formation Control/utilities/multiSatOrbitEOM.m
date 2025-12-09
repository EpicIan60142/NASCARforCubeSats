function dX = multiSatOrbitEOM(t, X, pConst, scConst, disturb)
% Function that integrates multiple satellite orbit EOMs at the same time
%   Inputs:
%       - t: Current integration time
%       - X: State vector of n stacked satellite states, like so:
%            [X_1; X_2; ...; X_n], where X_i = [X; Y; Z; Xdot; Ydot; Zdot]
%       - pConst: Planetary constants structure containing Ri, mu, and J2
%       - scConst: Vector of spacecraft parameter structures containing
%                  mass, coefficient of drag, and cross-sectional area
%       - disturb: Boolean of disturbances to include or ignore in the
%                  following order:
%           - J2: Include (true), ignore (false)
%           - Drag: Include (true), ignore (false)
%   Outputs:
%       - dX: Rate of change of the state vectors, like so:
%             [dX_1; dX_2; ...; dX_n]
%
%   By: Ian Faber, 08/26/2025
%

% Determine number of satellites in the combined state vector, each
% satellite has 6 states
nSats = length(X)/6;

dX = [];

% Calculate EOM for each satellite
for k = 1:nSats
    idx = (6*(k-1) + 1):(6*(k-1) + 6);
    X_i = X(idx);
    dX_i = orbitEOM(t, X_i, pConst, scConst(k), disturb);
    dX = [dX; dX_i];
end

end
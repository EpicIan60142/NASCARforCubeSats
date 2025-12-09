function dX = fullCWHEOM(t, X, const)
% Function that integrates the full Clohessy-Wiltshire-Hill equations 
% of motion between a chief spacecraft and a number of deputy spacecraft in
% the absence of perturbations.
%   Inputs:
%       - t: Current simulation time in seconds
%       - X: Formation state to integrate, organized as follows:
%            X = [X_chief; X_deputy1; X_deputy2; ... X_deputyN], where
%                 X_chief = [X; Y; Z; Xdot; Ydot; Zdot] in the inertial
%                           frame
%                 X_deputyN = [x; y; z; xDot; yDot; zDot] in the Hill frame
%       - const: Constants structure with a field for mu, the gravitational
%                parameter for the celestial body of interest
%   Outputs:
%       - dX: Rate of change vector for the formation as follows:
%             dX = [dX_chief; dX_deputy1; dX_deputy2; ... dX_deputyN]
%
%   By: Ian Faber, 09/16/2025
%

% Number of states being integrated
nStates = 6;

% Extract chief state
X_chief = X(1:nStates);

% Extract deputies and determine the number that exist
X_deputies = X(nStates+1:end);

nDeputies = length(X_deputies)/nStates;

% Calculate necessary intermediate quantities
    % Calculate DCM from inertial to Hill Frame
cart.rVec = X_chief(1:3);
cart.vVec = X_chief(4:6);
cart.mu = const.mu;
oe = convCart2ClassicOE(cart);

    % Orbit mean motion
n = sqrt(const.mu/(oe.a^3));

    % Chief orbit radius
r_c = norm(cart.rVec);

% Integrate chief state according to a simple point mass
rDDot_c = -(const.mu/(r_c^3))*cart.rVec;
dX_chief = [cart.vVec; rDDot_c];

% Integrate deputy states
dX_deputies = [];
for k = 1:nDeputies
        % Extract deputy information
    x = X_deputies((k-1)*nStates + 1);
    y = X_deputies((k-1)*nStates + 2);
    z = X_deputies((k-1)*nStates + 3);
    xDot = X_deputies((k-1)*nStates + 4);
    yDot = X_deputies((k-1)*nStates + 5);
    zDot = X_deputies((k-1)*nStates + 6);

    r_d = norm([r_c + x; y; z]);

        % Calculate accelerations
    xDDot = 2*n*yDot + 3*(n^2)*x; 
    yDDot = -2*n*xDot;
    zDDot = -(n^2)*z;

        % Append rate of change to accumulated vector
    dX_i = [xDot; yDot; zDot; xDDot; yDDot; zDDot];
    dX_deputies = [dX_deputies; dX_i];
end

% Assign output
dX = [dX_chief; dX_deputies];

end
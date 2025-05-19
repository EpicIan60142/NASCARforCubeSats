function cubesat = generateCubesat(uMax, axialThrust, startRing, name, marker, color)
% Function that generates CubeSats based on specified parameters
%   Inputs:
%       - uMax: Maximum thrust the CubeSat can produce
%       - axialThrust: Boolean indicating whether to split the specified
%                      max thrust up between 3 axial thrusters: 1 for x,
%                      one for y, and one for z axes
%       - startRing: Starting ring for the CubeSats to sample their initial
%                    state from
%       - name: What to name the CubeSat! Choose something fun ;)
%       - marker: What marker to plot the CubeSat as. Must be a valid
%                 'Marker' LineStyle option
%       - color: What color to plot the CubeSat trajectory as.
%   Outputs:
%       - cubesat: CubeSat output structure with the following fields:
%           - uMax: Max thrust, will either be a 3x1 vector or 1x1 scalar
%                   depending on the value of axialThrust:
%               - axialThrust == true: uMax = [uMax_x; uMax_y; uMax_z]
%               - axialThrust == false: uMax = uMax
%           - X0: Initial state for the current trajectory. At ring 1, 
%                 this is sampled from P0 like so:
%                 X0 = [X_sampled; Y_sampled; Z_sampled; 0; 0; 0]
%           - t0: Initial time for the current trajectory. At ring 1, 
%                 t0 = 0
%           - X: Trajectory of this Cubesat for the race course, organized
%                as follows:
%                X = [X, Y, Z, Xdot, Ydot, Zdot, px, py, pz, pxDot, pyDot,
%                     pzDot]
%           - t: Time at each point of this Cubesat's trajectory for the 
%                race course
%           - u: Control vector at each point of this Cubesat's trajectory
%                for the race course, organized as follows:
%                u = [u_x, u_y, u_z]
%           - initGuess: Initial guess for the parameters of this cube sat
%                        during each course segment, organized as follows:
%                        x0 = [p0; t_f]
%           - optParams: Parameters solved for this Cubesat during each
%                        course segment, organized as follows:
%                        optParams = [x_01, x_12, ..., x_lm1l], 
%                        where l is the number of rings in the course and 
%                        x = [p0; t_f]
%           - ringSeg: Stores what segment each set of parameters
%                      corresponds to, organized as follows:
%                      ringSeg = [ringSeg_1, ringSeg_2, ..., ringSeg_l], 
%                      where 
%                      ringSeg_i = [fromRingNum; toRingNum]
%           - tSeg: Stores the time that the cubesat started and finished
%                   each ring segment, organized as follows:
%                   tSeg = [tSeg_1, tSeg_2, ..., tSeg_l], where
%                   tSeg_i = [startRingSegTime; endRingSegTime]
%           - cost: Actual and minimum cost for the Cubesat to traverse
%                   each course segment, organized as follows:
%                   cost = [actualCost; minCost] 
%           - name: Specified cubesat name
%           - marker: Specified cubesat marker
%           - color: Specified cubesat color
%
%   By: Ian Faber, 04/19/2025
%
    % Define CubeSat Thrust
if axialThrust
    uMax_i = uMax/3;
    uMax = [uMax_i; uMax_i; uMax_i];
end

    % Define CubeSat starting position
        % Project starting ring into 3-D space
T = startRing.NR; % Get DCM for intertial coordinates from ring coordinates
T = T(:,1:2); % Get projection from nPrime and nDPrime vectors to XYZ space
P = T*startRing.S*T'; % 3sigma Covariance in 3D space

        % Sample starting position
sqrtP = zeros(size(P));
for k = 1:size(P,1)
    sqrtP(k,k) = sqrt(P(k,k));
end
X0_pos = startRing.center + sqrtP*randn(3,1)/3;
X0_vel = zeros(3,1);
X0 = [X0_pos; X0_vel];

    % Assign outputs
cubesat = struct("uMax", uMax, "X0", X0, "t0", 0, "X", [], "t", [], "u", [], ...
                 "initGuess", [], "optParams", [], "ringSeg", [], "tSeg", [], "cost", [], ...
                 "name", name, "marker", marker, "color", color);

end
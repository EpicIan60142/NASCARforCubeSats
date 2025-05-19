function ring = generateRing(a, b, dTheta, dPhi, d, lastRing)
% Function that generates a race course ring given a set of sizing and
% inter-ring distance parameters
%   Inputs:
%       - a: Semimajor axis of the ring to generate
%       - b: Semiminor axis of the ring to generate
%       - dTheta: Relative azimuth angle between the last ring's center and
%                 the new one's
%       - dPhi: Relative elevation angle between the last ring's center 
%               and the new one's
%       - d: Straight line distance from the center of the last ring to the
%            new one's
%       - lastRing: Ring structure for the last generated ring
%   Outputs:
%       - ring: Ring output structure with the following fields:
%           - center: Hill Frame coordinates of the center of the generated
%                     ring, organized as follows:
%                     [x_0; y_0; z_0]
%           - normal: Normal unit vector of the generated ring in the Hill
%                     Frame, organized as follows:
%                     n = [x_f - x_0; y_f - y_0; z_f - z_0], 
%                     normal = n/norm(n)
%           - S: Diagonal matrix encoding the shape of the generated ring
%                like so:
%                S = diag([a, b]);
%           - NR: DCM from ring frame to inertial frame, constructed as
%                 follows:
%                 NR = [normalPrime, normalDPrime, normal]
%           - params: Parameter structure with parameters that generated
%                     this ring, with the following fields from the 
%                     function input:
%               - a
%               - b
%               - theta = lastTheta + dTheta (converted to degrees)
%               - phi = lastPhi + dPhi (converted to degrees)
%               - d
%               - normalPrime
%               - normalDPrime
%               - lastRing
%
%   By: Ian Faber, 03/22/2025
%

%% Create ellipse size matrix
ring.S = diag([a, b]);

%% Generate new ring normal vector and center
    % Normal vector
thetaNew = deg2rad(lastRing.params.theta) + dTheta;
% thetaNew = theta;
phiNew = deg2rad(lastRing.params.phi) + dPhi;
% phiNew = phi;
n_x = abs(d)*cos(phiNew)*cos(thetaNew);
n_y = abs(d)*cos(phiNew)*sin(thetaNew);
n_z = abs(d)*sin(phiNew);
ring.normal = [n_x; n_y; n_z]/norm([n_x; n_y; n_z]);

normalPrime = [
                        sin(phiNew)*cos(thetaNew);
                        sin(phiNew)*sin(thetaNew);
                        -cos(phiNew)
                   ];

normalDPrime = cross(ring.normal, normalPrime);

    % Ring center
ring.center = lastRing.center + d*lastRing.normal;

    % DCM from ring frame to inertial frame
ring.NR = [normalPrime, normalDPrime, ring.normal];

%% Assign params structure
ring.params = struct("a", a, "b", b, "theta", rad2deg(thetaNew), "phi", rad2deg(phiNew), ...
                     "d", d, "normalPrime", normalPrime, "normalDPrime", normalDPrime, "lastRing", lastRing);


end





function ax = plotRing(ring, lineStyle)
% Function that plots a ring with a given center, sizing matrix and normal
% vector
%   Inputs:
%       - ring: Ring structure as defined in generateRing.m
%       - ringIdx: The number of the ring in the race course sequence, used
%                  for coloring
%       - lineStyle: Linestyle to plot the ring in. If empty, defaults to
%                    'k-'
%   Outputs:
%       - ax: Axis graphics handle for the ring
%
%   By: Ian Faber, 03/23/2025
%
%% Parse ring params
center = ring.center;
S = ring.S;

%% Extract linestyle
if isempty(lineStyle)
    lineStyle = 'k-';
end

%% Create theta vector for ellipse and pull out eigenvalues and eigenvectors
thetaCirc = linspace(0,2*pi,100);

[eigVec, Lambda] = eig(S);

%% Make template circle and scale it
circle = [cos(thetaCirc); sin(thetaCirc)];

circleScaled = (eigVec*sqrt(Lambda))*circle;

    % Add a z component at all zeroes
circleScaled = [circleScaled; zeros(size(thetaCirc))];

%% Rotate the z axis to coincide with the normal vector
%     % Define DCM from circle frame to normal frame, where the normal vector
%     % is the new z axis
% theta = deg2rad(ring.params.theta);
% phi = deg2rad(ring.params.phi);
% 
% NC = ...
%     [
%         sin(phi)*cos(theta), sin(phi)*sin(theta), -cos(phi); % normalPrime
%         -sin(theta),         cos(theta),          0          % normalDoublePrime
%         cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi);  % normal
%     ];

    % Apply transform
% circleRotated = NC'*circleScaled;
circleRotated = ring.NR*circleScaled;

%% Move center to coincide with the ring
circleTranslated = circleRotated + center;

%% Plot ring
ax = plot3(circleTranslated(1,:), circleTranslated(2,:), circleTranslated(3,:), lineStyle);

end
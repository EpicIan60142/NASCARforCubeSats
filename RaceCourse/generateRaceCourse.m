function rings = generateRaceCourse(courseParams)
% Function that generates a CubeSat race course given a number of specified
% course rules, i.e. sizing of rings, spacing and orientation between
% rings, and the number of rings in the course
%   Inputs:
%       - courseParams: Race course parameters structure with the following
%                       fields:
%           - semiMaj: Range of semi major axis sizes to choose from for
%                      each ring, specified as [min, max]
%           - semiMin: Range of semi minor axis sizes to choose from for
%                      each ring, specified as [min, max]
%           - dist: Range of distances between ring centers to choose from,
%                   specified as [min, max]
%           - azAng: Range of azimuth angles between ring centers to choose
%                    from in radians, specified as [min, max]
%           - elAng: Range of elevation angles between ring centers to
%                    choose from in radians, specified as [min, max]
%           - numRings: Range of the number of rings that make up the race
%                       course, specified as [min, max]
%   Outputs:
%       - rings: Race course ring structure array with a random number of
%                ring structures sampled from courseParams.numRings,
%                organized as follows:
%                [ring_1; ring_2; ...; ring_l], where l is the chosen
%                                               number of rings
%
%   By: Ian Faber, 03/22/2025
%

%% Parse courseParams
semiMaj = courseParams.semiMaj;
semiMin = courseParams.semiMin;
dist = courseParams.dist;
azAng = courseParams.azAng;
elAng = courseParams.elAng;
numRings = courseParams.numRings;

%% Choose number of rings to make
l = randi(numRings, 1);

%% Set up rings output
rings = [];

%% Create rings
for k = 1:l
    % Sample random parameters
        % Semimajor and semiminor axis
    a = min(semiMaj) + (max(semiMaj) - min(semiMaj))*rand(1);
    b = min(semiMin) + (max(semiMin) - min(semiMin))*rand(1);

        % Inter-ring distance
    d = min(dist) + (max(dist) - min(dist))*rand(1);

        % Azimuth and elevation angles
    if k <= l/2 % First half of the rings should turn left and down
        dTheta = min(azAng) + (0*max(azAng) - min(azAng))*rand(1);
        dPhi = min(elAng) + (0*max(elAng) - min(elAng))*rand(1);
    else % Second half of the rings should turn right and up
        dTheta = 0*min(azAng) + (max(azAng) - 0*min(azAng))*rand(1);
        dPhi = 0*min(elAng) + (max(elAng) - 0*min(elAng))*rand(1);
    end

        % Make rings
    if k == 1 % initial ring
            % Initial ring center
        center = [0; min(dist); 0];

            % Initial ring size
        a = max(semiMaj);
        b = max(semiMin);
        S = diag([a, b]);

            % Initial ring angles
        theta = pi/2;
        phi = 0;

            % Calculate initial normal vector and make DCM
        normal = [cos(phi)*cos(theta); cos(phi)*sin(theta); sin(phi)];
        normalPrime = [sin(phi)*cos(theta); sin(phi)*sin(theta); -cos(phi)];
        normalDPrime = cross(normal, normalPrime);

        NR = [normalPrime, normalDPrime, normal];

            % Create ring
        init = struct("center", center, "normal", normal, "S", S, "NR", NR); % Start first ring at the smallest allowed distance along the +y axis as the largest allowed ellipse
        init.params = struct("a", a, "b", b, "theta", rad2deg(theta), "phi", rad2deg(phi), ...
                             "d", 0, "normalPrime", normalPrime, "normalDPrime", normalDPrime, "lastRing", init);

        rings = [rings; init];
    else % All other rings
        rings = [rings; generateRing(a, b, dTheta, dPhi, d, rings(k-1))];
    end

end

end







%% ASEN 6014 Priject 1 Main Script
% By: Ian Faber, 10/07/2025

%% Housekeeping
clc; clear; close all;

addpath("Plotting/")
addpath("RaceCourse/")
addpath("Dynamics/")
addpath("../utilities/")

videoFolder = "Videos/";

%% Setup
    % Planetary constants structure
        % Point mass
pConst.mu = 398600.4415; % km^3/m^2
        % J2
pConst.J2 = 1.08264e-3;
pConst.Ri = 6378; % km
        % Drag
pConst.rho0 = 3.614e-4; % kg/km^3
pConst.r0 = 700 + pConst.Ri; % km
pConst.H = 88.667; % km

    % Race course settings
        % Ring geometry
semiMaj = [1, 5]/1000; % m->km, range of ring semi-major axes
semiMin = [1, 5]/1000; % m->km, range of ring semi-minor axes

        % Course geometry
interRingDist = [300, 350]/1000; % m, range of ring center separation distances
azimuthAngle = deg2rad([-90, 90]); % deg -> rad, relative azimuth angle range between rings
elevationAngle = deg2rad([-90, 90]); % deg -> rad, relative elevation angle range between rings
numRings = 13; % Num rings - 2, need to add initial and final rings

        % Course parameters structure
courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                      "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings);

        % Recover number of rings
numRings = numRings + 2;

        % Course origin orbit settings
T = 120*60; % 120 minute orbit period
oe_c.mu = pConst.mu;
oe_c.a = (pConst.mu*(T/(2*pi))^2)^(1/3);
oe_c.a = 8000;
oe_c.e = 0.18; % Decently eccentric course origin orbit
oe_c.i = deg2rad(45);
oe_c.RAAN = deg2rad(15);
oe_c.argPeri = deg2rad(-23);
oe_c.f = deg2rad(10);

    % ode45 settings
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

%% Task 1: Generate random race course
    % Determine if a random or seeded course is made
randomCourse = false;

if ~randomCourse
    seed = 3;
    rng(seed);
else
    rng("shuffle")
end

    % Generate race course
        % Intermediate rings
rings = generateRaceCourse(courseParams);

        % Starting ring
startRing = generateRing(5*max(semiMaj), 5*max(semiMin), 0, 0, -max(interRingDist), rings(1));
startRing.m = 100;
rings(1).params.lastRing = startRing;
rings = [startRing; rings];

        % End ring
endRing = generateRing(min(semiMaj), min(semiMin), 0, 0, min(interRingDist), rings(end));
endRing.m = 1;
rings = [rings; endRing];

    % Calculate race course origin and move center to be 0,0
centers = [];
for k = 1:length(rings)
    centers = [centers, rings.center];
end

newOrigin = mean(centers, 2);

for k = 1:length(rings)
    rings(k).center = rings(k).center - newOrigin;
    rings(k).params.lastRing.center = rings(k).params.lastRing.center - newOrigin;
end

    % Plot race course 
figNum = 1; 
titleText = sprintf("Generated Race Course");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
courseFig = plotCourse(rings, [], figNum, titleText, xLabel, yLabel, zLabel, [], []);

for k = 0:180
    view(-30 + 2*k, 25);
    drawnow;
end

    % Plot example ring

figNum = 2;
titleText = sprintf("Example Ring");
xLabel = sprintf("Radial [m]"); yLabel = sprintf("Along-Track [m]"); zLabel = sprintf("Cross-Track [m]");
exRingFig = plotExampleRing(rings, figNum, titleText, xLabel, yLabel, zLabel);
for k = 0:180
    view(-30 + 2*k, 35);
    drawnow;
end


%% Task 3/4: Develop integration scheme for all rings with J2 and drag and propagate over an average race
    % Compile initial ring states together
        % Inertial ring coordinates
cart_c = convClassicOE2Cart(oe_c);
X_c = [cart_c.rVec; cart_c.vVec];
X0_Inertial = X_c;
for k = 1:length(rings)
    ringRho = [rings(k).center; zeros(3,1)]; % All rings start at their center positions at rest
    ringX0 = convDeputyH2N(X_c, ringRho, pConst);
    X0_Inertial = [X0_Inertial; ringX0];
end

        % Hill ring coordinates
X0_Hill = X_c;
for k = 1:length(rings)
    X0_Hill = [X0_Hill; [rings(k).center; zeros(3,1)]];
end

    % Define integration time
tspan = 0:1:17*60; % Average race of 17 minutes

    % Loop over all disturbance cases: No disturbance, just J2, just drag,
    % J2 + drag
for k = 1:4
        % Determine disturbances for each ring
    switch k
        case 1
            J2_enable = false;
            drag_enable = false;
            titleText = sprintf("Relative Ring Trajectories with No Perturbation");
            movieTitle = videoFolder + sprintf("ASEN6014_Project1_NoPerturb.mp4");
            perturb = "No Perturbation";
        case 2
            J2_enable = true;
            drag_enable = false;
            titleText = sprintf("Relative Ring Trajectories with J_2 Perturbation");
            movieTitle = videoFolder + sprintf("ASEN6014_Project1_J2.mp4");
            perturb = "J_2 Perturbation";
        case 3
            J2_enable = false;
            drag_enable = true;
            titleText = sprintf("Relative Ring Trajectories with Drag Perturbation");
            movieTitle = videoFolder + sprintf("ASEN6014_Project1_Drag.mp4");
            perturb = "Drag Perturbation";
        case 4
            J2_enable = true;
            drag_enable = true;
            titleText = sprintf("Relative Ring Trajectories with J_2 + Drag Perturbation");
            movieTitle = videoFolder + sprintf("ASEN6014_Project1_J2Drag.mp4");
            perturb = "J_2 + Drag Perturbation";
    end

    disturb = [false; false]; % Course origin isn't ever subject to J2 or drag
    for kk = 1:length(rings)
        disturb = [disturb, [J2_enable; drag_enable]];
    end
    
        % Integrate ring motion forward in time in inertial coordinates
    [t_Inertial, X_Inertial] = ode45(@(t,X)multiSatOrbitEOMInertial(t,X,pConst,[rings(1);rings],disturb), tspan, X0_Inertial, opt);
    
        % Assign ring trajectories in the Hill frame
    rings_Inertial = rings;
    for kk = 1:length(rings_Inertial)
        ringX = [];
        for idx = 1:length(t_Inertial)
            ringX = [ringX, convDeputyN2H(X_Inertial(idx,1:6)', X_Inertial(idx,(kk*6 + 1):(kk*6 + 6))',pConst)];
        end
        rings_Inertial(kk).X = ringX;
        rings_Inertial(kk).t = t_Inertial;
    end
    
        % Integrate ring motion forward in time in Hill frame
    [t_Hill, X_Hill] = ode45(@(t,X)multiSatOrbitEOMRelative(t,X,pConst,[rings(1);rings],disturb), tspan, X0_Hill, opt);
    
        % Assign ring trajectories in the Hill frame
    rings_Hill = rings;
    for kk = 1:length(rings_Hill)
        rings_Hill(kk).X = X_Hill(:, (6*kk + 1):(6*kk + 6))';
        rings_Hill(kk).t = t_Hill;
    end
    
        % Plot numerically integrated rings
    figNum = (k-1)*3 + 3; 
    xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
    trajStyle = ["b-"; "r--"]; trajLabel = [sprintf("Ring Trajectories - Inertial mapped to Hill"); sprintf("Ring Trajectories - Integrated in Hill")];
    courseFig2 = plotCourse(rings_Inertial, rings_Hill, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);
    
        % Check answer
    figure((k-1)*3 + 4); tl = tiledlayout(3,1);
    title(tl, "Difference between Inertial and Hill EOM for each Ring")
    for kk = 1:length(rings_Inertial)
        nexttile(1);
            hold on; grid on;
            plot(t_Inertial, rings_Inertial(kk).X(1,:) - rings_Hill(kk).X(1,:))
            xlabel("Time [sec]"); ylabel("\Deltax");
        nexttile(2);
            hold on; grid on;
            plot(t_Inertial, rings_Inertial(kk).X(2,:) - rings_Hill(kk).X(2,:))
            xlabel("Time [sec]"); ylabel("\Deltay");
        nexttile(3);
            hold on; grid on;
            plot(t_Inertial, rings_Inertial(kk).X(3,:) - rings_Hill(kk).X(3,:))
            xlabel("Time [sec]"); ylabel("\Deltaz");
    end

    drawnow;

        % Animate ring motion
    animateRings(rings_Hill, true, movieTitle, (k-1)*3 + 5, perturb);
end




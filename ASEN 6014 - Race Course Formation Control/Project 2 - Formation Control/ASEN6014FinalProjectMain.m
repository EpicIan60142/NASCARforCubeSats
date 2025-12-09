%% ASEN 6014 Final Project Main Script
% By: Ian Faber, 11/19/2025

%% Housekeeping
clc; clear; close all;

addpath("Controls\")
addpath("Dynamics\")
addpath("Plotting\")
addpath("RaceCourse\")
addpath("../utilities/")

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
T = 180*60; % 180 minute orbit period
oe_c.mu = pConst.mu;
oe_c.a = (pConst.mu*(T/(2*pi))^2)^(1/3);
oe_c.e = 0.05; % Decently eccentric course origin orbit
oe_c.i = deg2rad(45);
oe_c.RAAN = deg2rad(15);
oe_c.argPeri = deg2rad(-23);
oe_c.f = deg2rad(10);

    % ode45 settings
opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);

    % Video saving
videoFolder = "Videos/";

%% Task 0: Generate a desired race course formation
    % Determine if a random course is made
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
endRing = generateRing(5*max(semiMaj), 5*max(semiMin), 0, 0, min(interRingDist), rings(end));
endRing.m = 100;
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

titleText = sprintf("Generated Race Course");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
plotCourse(rings, 420, titleText, xLabel, yLabel, zLabel, [], []);

return

%% Task 1a: Design and implement a race course formation feedback control law
    % Define gains, goal is to go from lead-follower to end position within
    % 10 minutes
kConst_course.K1 = 2e-5;
kConst_course.K2 = 5e-3;

    % Get initial state of the race course origin
cart_c = convClassicOE2Cart(oe_c);
X_c0 = [cart_c.rVec; cart_c.vVec];

    % Test ring formation control on a random ring from the race course
X_d0 = convDeputyH2N(X_c0, zeros(6,1), pConst); % Start at the course origin
X0_test = [X_c0; X_d0];
tspan_test = 0:10:T; % Ideally want rings to reach their positions within an hour
ringIdx = randi(length(rings),1);

[t_test, X_test] = ode23t(@(t,X)ringFormationFeedbackControl(t, X, rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), tspan_test, X0_test);
[~, u_test] = cellfun(@(t,X)ringFormationFeedbackControl(t, X.', rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), num2cell(t_test), num2cell(X_test,2), 'uni', 0);
u_test = cellfun(@(x)x', u_test, 'UniformOutput', false);
u_test = cell2mat(u_test);

    % Convert deputy from Inertial frame to Hill frame
X_test_Hill = [];
for k = 1:length(t_test)
    X_test_Hill = [X_test_Hill; convDeputyN2H(X_test(k, 1:6)', X_test(k, 7:12)', pConst)'];
end

    % Assign trajectory arrays
rings(ringIdx).X = X_test_Hill;
rings(ringIdx).t = t_test;
rings(ringIdx).u = u_test;
rings(ringIdx).X_r = repmat([rings(ringIdx).center; zeros(3,1)]', length(t_test), 1);
rings(ringIdx).oe = NaN(length(t_test), 6); % Not controlling orbit elements
rings(ringIdx).oe_r = NaN(length(t_test), 6);

%% Task 1b: Design and implement a lead-follower formation feedback control law
    % Define lead follower characteristics
doe_r.da = 0; % km
doe_r.de = 0;
doe_r.di = deg2rad(0); % deg -> rad
doe_r.dRAAN = deg2rad(0); % deg -> rad
doe_r.dargPeri = deg2rad(0); % deg -> rad
doe_r.dM0 = deg2rad(2e-4); % deg -> rad

kConst_lead.P11 = 4e-3;
kConst_lead.P22 = 8e-3;
kConst_lead.P33 = 4e-3;

kConst_lead.K11 = 8e-3;
kConst_lead.K22 = 4e-3;
kConst_lead.K33 = 8e-3;
kConst_lead.K44 = 8e-3;
kConst_lead.K55 = 4e-3;
kConst_lead.K66 = 4e-3;

    % Add another orbit of time
tspan_test = (0:10:3*T) + t_test(end); % Want rings to stow within two orbits

    % Update initial state
X0_test = X_test(end,:)';

    % Run lead-follower control
[t_test, X_test] = ode45(@(t,X)leadFollowerFeedbackControl(t, X, doe_r, kConst_lead, pConst, rings(ringIdx)), tspan_test, X0_test, opt);
[~, u_test, doe_test, oe_d_test, oe_c_test, oe_r_test] = cellfun(@(t,X)leadFollowerFeedbackControl(t, X.', doe_r, kConst_lead, pConst, rings(ringIdx)), num2cell(t_test), num2cell(X_test,2), 'uni', 0);
u_test = cellfun(@(x)x', u_test, 'UniformOutput', false);
u_test = cell2mat(u_test);
doe_test = cellfun(@(x)x',doe_test,'uni',0);
doe_test = cell2mat(doe_test);
oe_d_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_d_test, 'uni', 0);
oe_d_test = cell2mat(oe_d_test);
oe_c_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_c_test, 'uni', 0);
oe_c_test = cell2mat(oe_c_test);
oe_r_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_r_test, 'uni', 0);
oe_r_test = cell2mat(oe_r_test);

    % Convert deputy from Inertial frame to Hill frame
X_test_Hill = [];
for k = 1:length(t_test)
    X_test_Hill = [X_test_Hill; convDeputyN2H(X_test(k, 1:6)', X_test(k, 7:12)', pConst)'];
end

    % Assign trajectory arrays
rings(ringIdx).X = [rings(ringIdx).X; X_test_Hill];
rings(ringIdx).t = [rings(ringIdx).t; t_test];
rings(ringIdx).u = [rings(ringIdx).u; u_test];
rings(ringIdx).X_r = [rings(ringIdx).X_r; NaN(length(t_test), 6)]; % Not controlling Hill coords
rings(ringIdx).oe = [rings(ringIdx).oe; oe_d_test];
rings(ringIdx).oe_r = [rings(ringIdx).oe_r; oe_r_test];

%% Task 1c: Plot results of both controllers
    % Plot results
figNum = 1;
titleText = sprintf("Test of ring formation control");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
trajStyle = "b-"; trajLabel = sprintf("Ring Trajectory");
plotCourse(rings, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);

    % Plot control effort
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Test Control Effort vs. Time");
nt = nexttile; ax = [ax; nt];
    semilogy(rings(ringIdx).t, abs(rings(ringIdx).u(:,1)), 'b-');
    hold on; grid on;
    xlabel("Time [sec]"); ylabel("u_r [km/s^2]");
nt = nexttile; ax = [ax; nt];  
    semilogy(rings(ringIdx).t, abs(rings(ringIdx).u(:,2)), 'b-');
    hold on; grid on;
    xlabel("Time [sec]"); ylabel("u_\theta [km/s^2]");
nt = nexttile; ax = [ax; nt];
    semilogy(rings(ringIdx).t, abs(rings(ringIdx).u(:,3)), 'b-');
    hold on; grid on;
    xlabel("Time [sec]"); ylabel("u_h [km/s^2]");
linkaxes(ax, 'x');

    % Plot ring trajectory
figure; tl = tiledlayout(3,2);
title(tl, "Test Ring Trajectory vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,1), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,1), 'r--');
    xlabel("Time [sec]"); ylabel("x [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,4), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,4), 'r--');
    xlabel("Time [sec]"); ylabel("xDot [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,2), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,2), 'r--');
    xlabel("Time [sec]"); ylabel("y [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,5), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,5), 'r--');
    xlabel("Time [sec]"); ylabel("yDot [km/s]");
    legend("Trajectory", "Reference", 'location', 'eastoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,3), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,3), 'r--');
    xlabel("Time [sec]"); ylabel("z [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).X(:,6), 'b-');
    plot(rings(ringIdx).t, rings(ringIdx).X_r(:,6), 'r--');
    xlabel("Time [sec]"); ylabel("zDot [km/s]");
linkaxes(ax, 'x');

    % Plot ring elements
figure; tl = tiledlayout(3,2);
title(tl, "Test Ring Orbit Element Differences vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).oe(:,1) - rings(ringIdx).oe_r(:,1), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,1), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta a_d [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rings(ringIdx).oe(:,2) - rings(ringIdx).oe_r(:,2), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,2), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta e_d");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,3) - rings(ringIdx).oe_r(:,3)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,3), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta i_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,4) - rings(ringIdx).oe_r(:,4)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,4), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta \Omega_d [deg]");
    legend("Trajectory", "Reference", 'location', 'eastoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,5) - rings(ringIdx).oe_r(:,5)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,5), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta \omega_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(rings(ringIdx).t, rad2deg(rings(ringIdx).oe(:,6) - rings(ringIdx).oe_r(:,6)), 'b-');
    plot(rings(ringIdx).t, 0*rings(ringIdx).oe_r(:,6), 'r--');
    xlabel("Time [sec]"); ylabel("\Delta M_d [deg]");
linkaxes(ax, 'x');

drawnow;

%     % Animate ring
% animateRings(rings, true, videoFolder + "ControllerTest.mp4", 5, "");

%% Task 2: Design and implement a deployment sequence from lead-follower to race course formation
    % Clear reused parameters
clear doe_r oe_d X_d0

    % Clear the example ring trajectory
rings(ringIdx).t = [];
rings(ringIdx).X = [];
rings(ringIdx).u = [];
rings(ringIdx).X_r = [];
rings(ringIdx).oe = [];
rings(ringIdx).oe_r = [];

    % Sort the intermediate rings by increasing goal distance from the origin
distances = [];
for k = 2:length(rings)-1
    distances = [distances; [norm(rings(k).center), k]];
end

[sortDist, sortIdx] = sort(distances);
sortDist(:,2) = distances(sortIdx(:,1), 2);
sortIdx(:,2) = [];

    % Assign initial mean anomalies
minDM0 = 9e9;
maxDM0 = 0;
numPos = 0; % Number of rings starting with a positive dM0
numNeg = 0; % Number of rings starting with a negative dM0
deployOrder = [];
for k = 1:length(sortDist)
        % Define desired orbit element differences for this ring
    doe_r.da = 0;
    doe_r.de = 0;
    doe_r.di = 0;
    doe_r.dRAAN = 0;
    doe_r.dargPeri = 0;
    if mod(k,2) % Fan out rings to each side of the origin at the start
        doe_r.dM0 = deg2rad(k*1e-4 + 0.5e-4);
        numPos = numPos + 1;
    else
        doe_r.dM0 = deg2rad(k*-1e-4 + 0.5e-4);
        numNeg = numNeg + 1;
    end

    if doe_r.dM0 > maxDM0
        maxDM0 = doe_r.dM0;
    end

    if doe_r.dM0 < minDM0
        minDM0 = doe_r.dM0;
    end

    oe_d = oe_c;
    M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
    M_r = M_c + doe_r.dM0;
    oe_d.f = convE2f(convM2E(M_r, oe_d.e, false), oe_d.e);
    
        % Assign resulting cartesian state as the initial inertial state of
        % the ring
    cart_d = convClassicOE2Cart(oe_d);
    % rings(sortDist(k,2)).X0 = convDeputyN2H(X_c0, [cart_d.rVec; cart_d.vVec], pConst);
    rings(sortDist(k,2)).X0 = [cart_d.rVec; cart_d.vVec];
    % rings(sortDist(k,2)).X = convDeputyN2H(X_c0, rings(sortDist(k,2)).X0, pConst)';
    % rings(sortDist(k,2)).t = 0;

        % Assign initial lead follower desired elements
    rings(sortDist(k,2)).doe_r = doe_r;

        % Assign deployment order
    deployOrder = [deployOrder; sortDist(k,2)];
end

    % Assign start and end ring starting positions at the extremes of the
    % lead-follower configuration
        % Start ring - position on negative along-track side
doe_r.dM0 = minDM0 - 2*deg2rad(2e-4);
numNeg = numNeg + 1;
oe_d = oe_c;
M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
M_r = M_c + doe_r.dM0;
oe_d.f = convE2f(convM2E(M_r, oe_d.e, false), oe_d.e);
cart_d = convClassicOE2Cart(oe_d);
% rings(1).X0 = convDeputyN2H(X_c0, [cart_d.rVec; cart_d.vVec], pConst);
rings(1).X0 = [cart_d.rVec; cart_d.vVec];
% rings(1).X = convDeputyN2H(X_c0, rings(1).X0, pConst)';
% rings(1).t = 0;
rings(1).doe_r = doe_r;
deployOrder = [deployOrder; 1];

        % End ring - position on positive along-track axis
doe_r.dM0 = maxDM0 + 2*deg2rad(2e-4);
numPos = numPos + 1;
oe_d = oe_c;
M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
M_r = M_c + doe_r.dM0;
oe_d.f = convE2f(convM2E(M_r, oe_d.e, false), oe_d.e);
cart_d = convClassicOE2Cart(oe_d);
% rings(end).X0 = convDeputyN2H(X_c0, [cart_d.rVec; cart_d.vVec], pConst);
rings(end).X0 = [cart_d.rVec; cart_d.vVec];
% rings(end).X = convDeputyN2H(X_c0, rings(end).X0, pConst)';
% rings(end).t = 0;
rings(end).doe_r = doe_r;
deployOrder = [deployOrder; length(rings)];

    % Flip order to deploy the last added first
deployOrder = flip(deployOrder);

    % Command the rings to deploy
fprintf("Deploying race course...\n");
deployDelay = 500;
for k = 1:length(deployOrder)
        % Get the ring to deploy
    ringIdx = deployOrder(k);

    fprintf("\tDeploying ring %.0f! (%.0f/%.0f)\n",deployOrder(k),k,length(deployOrder));

        % Propagate controller
            % Result reporting time interval
    dt = 10;
            % Define deployment time for this ring
    tspan = (1*T + k*deployDelay):dt:3*T;
            % Start with orbit element control until deployment time
    if tspan(1) ~= 0 
        tspan_elem = 0:dt:tspan(1);
        cart_c = convClassicOE2Cart(oe_c);
        X_c0 = [cart_c.rVec; cart_c.vVec];
        [t_elem, X_elem] = ode45(@(t,X)leadFollowerFeedbackControl(t, X, rings(ringIdx).doe_r, kConst_lead, pConst, rings(ringIdx)), tspan_elem, [X_c0; rings(ringIdx).X0], opt);
        [~, u_elem, doe_elem, oe_d_elem, oe_c_elem, oe_r_elem] = cellfun(@(t,X)leadFollowerFeedbackControl(t, X.', rings(ringIdx).doe_r, kConst_lead, pConst, rings(ringIdx)), num2cell(t_elem), num2cell(X_elem,2), 'uni', 0);
        u_elem = cellfun(@(x)x', u_elem, 'UniformOutput', false);
        u_elem = cell2mat(u_elem);
        doe_elem = cellfun(@(x)x',doe_elem,'uni',0);
        doe_elem = cell2mat(doe_elem);
        oe_d_elem = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_d_elem, 'uni', 0);
        oe_d_elem = cell2mat(oe_d_elem);
        % oe_c_test = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_c_test, 'uni', 0);
        % oe_c_test = cell2mat(oe_c_test);
        oe_r_elem = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_r_elem, 'uni', 0);
        oe_r_elem = cell2mat(oe_r_elem);

            % Convert deputy from Inertial frame to Hill frame
        X_elem_Hill = [];
        for kk = 1:length(t_elem)
            X_elem_Hill = [X_elem_Hill; convDeputyN2H(X_elem(kk, 1:6)', X_elem(kk, 7:12)', pConst)'];
        end
        
            % Assign trajectory arrays
        rings(ringIdx).X = [rings(ringIdx).X; X_elem_Hill];
        rings(ringIdx).t = [rings(ringIdx).t; t_elem];
        rings(ringIdx).u = [rings(ringIdx).u; u_elem];
        rings(ringIdx).X_r = [rings(ringIdx).X_r; NaN(length(t_elem), 6)]; % Not controlling Hill coords
        rings(ringIdx).oe = [rings(ringIdx).oe; oe_d_elem];
        rings(ringIdx).oe_r = [rings(ringIdx).oe_r; oe_r_elem];
    end

            % Run formation control after deployment time starts
    if tspan(1) ~= 0
        oe_c0 = oe_c_elem{end};
        oe_c0.mu = pConst.mu;
        cart_c0 = convClassicOE2Cart(oe_c0);
        X_c0 = [cart_c0.rVec; cart_c0.vVec];
        rings(ringIdx).X0 = convDeputyH2N(X_c0, rings(ringIdx).X(end,:)', pConst);
    end
    [t, X] = ode23t(@(t,X)ringFormationFeedbackControl(t, X, rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), tspan, [X_c0; rings(ringIdx).X0], opt);
    [~, u] = cellfun(@(t,X)ringFormationFeedbackControl(t, X.', rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), num2cell(t), num2cell(X,2), 'uni', 0);
    u = cellfun(@(x)x', u, 'UniformOutput', false);
    u = cell2mat(u);
    
        % Convert deputy from Inertial frame to Hill frame
    X_Hill = [];
    for k = 1:length(t)
        X_Hill = [X_Hill; convDeputyN2H(X(k, 1:6)', X(k, 7:12)', pConst)'];
    end
    
        % Assign trajectory
    rings(ringIdx).X = [rings(ringIdx).X; X_Hill];
    rings(ringIdx).t = [rings(ringIdx).t; t];
    rings(ringIdx).u = [rings(ringIdx).u; u];
    rings(ringIdx).X_r = [rings(ringIdx).X_r; repmat([rings(ringIdx).center; zeros(3,1)]', length(t), 1)];
    rings(ringIdx).oe = [rings(ringIdx).oe; NaN(length(t), 6)]; % Not controlling orbit elements
    rings(ringIdx).oe_r = [rings(ringIdx).oe_r; NaN(length(t), 6)];
end

fprintf("\nRace course deployed, let's race!!!\n\n")

%     % Plot trajectories
% titleText = sprintf("Race Course Ring Deployment");
% xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
% trajStyle = "b-"; trajLabel = sprintf("Ring Trajectory");
% plotCourse(rings, 6, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);
% 
%     % Animate deployment
% animateRings(rings, true, videoFolder + "RingDeployment.mp4", 7, "");

%% Task 3: Design and implement a stowing sequence from the race course to lead-follower formation
    % Stowing order is in reverse of deployment order
stowOrder = flip(deployOrder);

    % Capture the initial state of the chief after formation control
X_c0_form = X(end,1:6)';

    % Command the rings to stow
fprintf("Stowing race course...\n");
stowDelay = 1000;
for k = 1:length(stowOrder)
        % Get the ring to deploy
    ringIdx = stowOrder(k);

    fprintf("\tStowing ring %.0f! (%.0f/%.0f)\n",stowOrder(k),k,length(stowOrder));

        % Propagate controller
            % Result reporting time interval
    dt = 10;
            % Define deployment time for this ring
    if k == 1
        oldTSpan = tspan;
    end

    tspan = (oldTSpan(end) + k*stowDelay):dt:(oldTSpan(end) + 7*T);
            % Start with race course formation control until stowing time
    if tspan(1) ~= oldTSpan(end)
        tspan_course = oldTSpan(end):dt:tspan(1);
        X_c0 = X_c0_form;
        rings(ringIdx).X0 = convDeputyH2N(X_c0, rings(ringIdx).X(end,:)', pConst);
        [t_course, X_course] = ode15s(@(t,X)ringFormationFeedbackControl(t, X, rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), tspan_course, [X_c0; rings(ringIdx).X0], opt);
        [~, u_course] = cellfun(@(t,X)ringFormationFeedbackControl(t, X.', rings(ringIdx).center, kConst_course, pConst, rings(ringIdx)), num2cell(t_course), num2cell(X_course,2), 'uni', 0);
        u_course = cellfun(@(x)x', u_course, 'UniformOutput', false);
        u_course = cell2mat(u_course);
        
            % Convert deputy from Inertial frame to Hill frame
        X_elem_Hill = [];
        for kk = 1:length(t_course)
            X_elem_Hill = [X_elem_Hill; convDeputyN2H(X_course(kk, 1:6)', X_course(kk, 7:12)', pConst)'];
        end
        
            % Assign trajectory arrays
        rings(ringIdx).X = [rings(ringIdx).X; X_elem_Hill];
        rings(ringIdx).t = [rings(ringIdx).t; t_course];
        rings(ringIdx).u = [rings(ringIdx).u; u_course];
        rings(ringIdx).X_r = [rings(ringIdx).X_r; repmat([rings(ringIdx).center; zeros(3,1)]', length(t_course), 1)]; % Not controlling Hill coords
        rings(ringIdx).oe = [rings(ringIdx).oe; NaN(length(t_course), 6)]; % Not controlling orbit elements
        rings(ringIdx).oe_r = [rings(ringIdx).oe_r; NaN(length(t_course), 6)];
    end

            % Run lead follower control after stow time starts
    if tspan(1) ~= oldTSpan(end)
        X_c0 = X_course(end,1:6)';
    else
        X_c0 = X_c0_form;
    end
    rings(ringIdx).X0 = convDeputyH2N(X_c0, rings(ringIdx).X(end,:)', pConst);
    [t, X] = ode45(@(t,X)leadFollowerFeedbackControl(t, X, rings(ringIdx).doe_r, kConst_lead, pConst, rings(ringIdx)), tspan, [X_c0; rings(ringIdx).X0], opt);
    [~, u, doe, oe_d, oe_c, oe_r] = cellfun(@(t,X)leadFollowerFeedbackControl(t, X.', rings(ringIdx).doe_r, kConst_lead, pConst, rings(ringIdx)), num2cell(t), num2cell(X,2), 'uni', 0);
    u = cellfun(@(x)x', u, 'UniformOutput', false);
    u = cell2mat(u);
    doe = cellfun(@(x)x',doe,'uni',0);
    doe = cell2mat(doe);
    oe_d = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_d, 'uni', 0);
    oe_d = cell2mat(oe_d);
    oe_r = cellfun(@(x)[x.a; x.e; x.i; x.RAAN; x.argPeri; convE2M(convf2E(x.f,x.e),x.e)]', oe_r, 'uni', 0);
    oe_r = cell2mat(oe_r);
    
        % Convert deputy from Inertial frame to Hill frame
    X_Hill = [];
    for k = 1:length(t)
        X_Hill = [X_Hill; convDeputyN2H(X(k, 1:6)', X(k, 7:12)', pConst)'];
    end
    
        % Assign trajectory
    rings(ringIdx).X = [rings(ringIdx).X; X_Hill];
    rings(ringIdx).t = [rings(ringIdx).t; t];
    rings(ringIdx).u = [rings(ringIdx).u; u];
    rings(ringIdx).X_r = [rings(ringIdx).X_r; NaN(length(t), 6)]; % Not controlling Hill coordinates
    rings(ringIdx).oe = [rings(ringIdx).oe; oe_d];
    rings(ringIdx).oe_r = [rings(ringIdx).oe_r; oe_r];
end

fprintf("\n\nRace Course stowed, come back soon!!!\n\n");

%% Task 4: Plot full deployment - sustain - stow sequence
    % Plot trajectories
titleText = sprintf("Race Course Ring Deployment and Stowing");
xLabel = sprintf("Radial [km]"); yLabel = sprintf("Along-Track [km]"); zLabel = sprintf("Cross-Track [km]");
trajStyle = "b-"; trajLabel = sprintf("Ring Trajectory");
plotCourse(rings, 8, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel);

    % Plot control effort
markerSize = 5; colorStyle = 'cool';
figure; tl = tiledlayout(3,1); ax = [];
title(tl, "Race Course Control Effort vs. Time");
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        a = scatter(rings(k).t, abs(rings(k).u(:,1)), markerSize, k*ones(size(rings(k).t)), 'filled');
    end
    set(gca, 'yscale', 'log'); colormap(colorStyle); 
    xlabel("Time [sec]"); ylabel("u_r [km/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        a = scatter(rings(k).t, abs(rings(k).u(:,2)), markerSize, k*ones(size(rings(k).t)), 'filled');
    end
    set(gca, 'yscale', 'log'); colormap(colorStyle);
    cbar = colorbar; cbar.Label.String = "Ring number"; cbar.Location = 'eastoutside';
    xlabel("Time [sec]"); ylabel("u_\theta [km/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        a = scatter(rings(k).t, abs(rings(k).u(:,3)), markerSize, k*ones(size(rings(k).t)), 'filled');
    end
    set(gca, 'yscale', 'log'); colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("u_h [km/s^2]");
linkaxes(ax, 'x');

    % Plot ring trajectory
figure; tl = tiledlayout(3,2);
title(tl, "Race Course Ring Trajectory vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).X(:,1), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, rings(k).X_r(:,1), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("x [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).X(:,4), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, rings(k).X_r(:,4), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("xDot [km/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).X(:,2), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, rings(k).X_r(:,2), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("y [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        ring = scatter(rings(k).t, rings(k).X(:,5), markerSize, k*ones(size(rings(k).t)), 'filled');
        ref = plot(rings(k).t, rings(k).X_r(:,5), 'r--');
    end
    cbar = colorbar; cbar.Label.String = "Ring number"; cbar.Location = 'eastoutside';
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("yDot [km/s]");
    legend([ring, ref], ["Trajectory", "Reference"], 'location', 'eastoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).X(:,3), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, rings(k).X_r(:,3), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("z [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on; grid minor;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).X(:,6), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, rings(k).X_r(:,6), 'r--');
    end
    xlabel("Time [sec]"); ylabel("zDot [km/s]");
linkaxes(ax, 'x');

    % Plot ring elements
figure; tl = tiledlayout(3,2);
title(tl, "Race Course Ring Orbit Element Differences vs. Time")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).oe(:,1) - rings(k).oe_r(:,1), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, 0*rings(k).oe_r(:,1), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("\Delta a_d [km]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).oe(:,2) - rings(k).oe_r(:,2), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, 0*rings(k).oe_r(:,2), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("\Delta e_d");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).oe(:,3) - rings(k).oe_r(:,3), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, 0*rings(k).oe_r(:,3), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("\Delta i_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(rings)
        ring = scatter(rings(k).t, rings(k).oe(:,4) - rings(k).oe_r(:,4), markerSize, k*ones(size(rings(k).t)), 'filled');
        ref = plot(rings(k).t, 0*rings(k).oe_r(:,4), 'r--');
    end
    cbar = colorbar; cbar.Label.String = "Ring number"; cbar.Location = 'eastoutside';
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("\Delta \Omega_d [deg]");
    legend([ring, ref], ["Trajectory", "Reference"], 'location', 'eastoutside')
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).oe(:,5) - rings(k).oe_r(:,5), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, 0*rings(k).oe_r(:,5), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("\Delta \omega_d [deg]")
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    for k = 1:length(rings)
        scatter(rings(k).t, rings(k).oe(:,6) - rings(k).oe_r(:,6), markerSize, k*ones(size(rings(k).t)), 'filled');
        plot(rings(k).t, 0*rings(k).oe_r(:,6), 'r--');
    end
    colormap(colorStyle);
    xlabel("Time [sec]"); ylabel("\Delta M_d [deg]");
linkaxes(ax, 'x');

    % Animate deployment
animateRings(rings, true, videoFolder + "RingDeploymentAndStowing.mp4", 12, "");

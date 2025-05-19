%% ASEN 6020 Final Project Main Script
% By: Ian Faber

%% Housekeeping
clc; clear; close all;

%% Logistics
addpath(genpath(".\CubeSats"));
addpath(genpath(".\RaceCourse"));
addpath(genpath(".\Plotting"));

% %% Enable parallel computing
% parpool;

%% Ask user for input
menuMsg = sprintf("\n\t\tWelcome to the CubeSat racer! What would you like to do?\n\n" + ...
                  "\t1. Generate and run a new race course (time intensive, ~30 min to 1 hour!)\n" + ...
                  "\t2. Analyze an old run\n" + ...
                  "\t3. Animate an old run\n" + ...
                  "\t4. Exit\n\n" + ...
                  "\tChoice (enter number of option to do): ");
choice = input(menuMsg);

%% Decide what to do
switch choice
    case 1
        %% Setup Problem Pameters
        menuMsg = sprintf("\n\tRun dev course or generate new course?\n" + ...
                          "\t1. Run dev course (different thrust configs)\n" + ...
                          "\t2. Generate new course (same thrust configs)\n" + ...
                          "\tChoice (enter number of option to do): ");
        choice = input(menuMsg);

        if choice == 1
            randomCourse = false;
            fprintf("Running dev race course, this will take some time!\n")
        else
            randomCourse = true;
            fprintf("Generating and solving new race course, this will take some time!\n")
        end

            % Define ring semi-major and minor axis ranges
        semiMaj = [1, 5]; % m
        semiMin = [1, 5]; % m
        
            % Define inter-ring parameter ranges
        interRingDist = [300, 350]; % m
        azimuthAngle = deg2rad([-90, 90]); % deg -> rad
        elevationAngle = deg2rad([-90, 90]); % deg -> rad
        
            % Define number of rings range
        numRings = [15, 25];
        
            % Define course origin orbit
        T = 180*60; % 180 minute circular orbit
        n = 2*pi/T;
        
            % Define course parameters structure
        courseParams = struct("semiMaj", semiMaj, "semiMin", semiMin, "dist", interRingDist, ...
                              "azAng", azimuthAngle, "elAng", elevationAngle, "numRings", numRings, ...
                              "n", n);
        
            % ODE45 options
        opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
        
        %% Generate race course
            % Set rng seed for race course generation
        if ~randomCourse
            seed = 2*69420; % cool options: 0, 2, 3, 4
            rng(seed);
        else
            rng("shuffle");
        end
        
            % Make course
        rings = generateRaceCourse(courseParams);
        
            % Define CubeSat starting ring as the maximum distance behind the start
            % ring at a 3 sigma covariance of 5x the largest semi-major/minor axis
        startRing = generateRing(50*max(semiMaj), 50*max(semiMin), 0, 0, -max(interRingDist), rings(1));
        rings(1).params.lastRing = startRing;

            % Define CubeSat ending ring as the minimum distance after the last 
            % ring at a covariance corresponding to the smallest possible ring
        endRing = generateRing(min(semiMaj),min(semiMin),0,0,min(interRingDist),rings(end));
        rings = [rings; endRing];
        
            % Move race course so origin is the mean of all ring centers
        centers = [];
        for k = 1:length(rings)
            centers = [centers, rings.center];
        end
        
        newOrigin = mean(centers, 2);
        
        startRing.center = startRing.center - newOrigin;
        startRing.params.lastRing.center = startRing.params.lastRing.center - newOrigin;
        for k = 1:length(rings)
            rings(k).center = rings(k).center - newOrigin;
            rings(k).params.lastRing.center = rings(k).params.lastRing.center - newOrigin;
        end
        endRing.center = endRing.center - newOrigin;
        endRing.params.lastRing.center = endRing.params.lastRing.center - newOrigin;
        
        %% Generate CubeSats
            % Max thrust based on
            % https://www.cubesatshop.com/product/cluster-ifm-nano-thruster-smallsats/
            % - 10-500 uN. Boosted by 1e3 for this problem to make things
            % more interesting!
        numSats = 4;
        names = ["Eeny", "Meeny", "Miny", "Mo"]; % Names :P Mo has attitude control system problems, and can only use x-y-z thrusters!
        markers = ["o", "^", "square", "diamond"]; % Markers for plotting
        colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]; % Colors for plotting
        if ~randomCourse
            thrusts = [100, 250, 500, 500]*(10^-3); % Maximum thrust values
            thrustConfigs = [false, false, false, true]; % Whether thrust is split amongst x-y-z thrusters or a general direction
        else
            thrusts = 500*10^-3*ones(1,4); % Make cubesats with max thrust
            thrustConfigs = false(1,4); % No thrust configuration shenanigans
        end

            % Make CubeSats
        cubesats = [];
        for k = 1:numSats
            cubesats = [cubesats; generateCubesat(thrusts(k), thrustConfigs(k), startRing, names(k), markers(k), colors(k))];
        end
        
        %% Show off course and cubesat starting positions
        figNum = 420;
        titleText = sprintf("Generated Race Course");
        xLabel = sprintf("Radial [m]"); yLabel = sprintf("Along Track [m]"); zLabel = sprintf("Cross Track [m]");
        courseFig = plotCourse(startRing, rings, endRing, cubesats, figNum, titleText, xLabel, yLabel, zLabel);

        for k = 0:180
            view(-30 + 2*k, 35);
            drawnow
        end
        
        %% Run the race course!
        debug = [false; false; false]; % iteration plotting+display; trajectory plotting; disable progress sequence
        for k = 1:length(cubesats)
            fprintf("\nCubesat %s is starting the race course!\n", cubesats(k).name)
            for kk = 1:length(rings) % Loop through each course segment
                    % Determine which cost function and constraints to
                    % apply
                if kk == length(rings) % Final problem, come to a stop in minimum time
                    isFinal = true;
                else % Initial/intermediate problem, aim for rings
                    isFinal = false;
                end

                    % Solve trajectory according to problem
                [t,X,u,optParams,initGuess,cost] = solveTrajectory(cubesats(k), rings(kk), courseParams, opt, isFinal, debug);
                cubesats(k).X = [cubesats(k).X; X];
                cubesats(k).t = [cubesats(k).t; t];
                cubesats(k).u = [cubesats(k).u; u];
                cubesats(k).optParams = [cubesats(k).optParams, optParams]; % Parameters solved by fmincon
                cubesats(k).initGuess = [cubesats(k).initGuess, initGuess]; % Initial guess for parameters
                cubesats(k).ringSeg = [cubesats(k).ringSeg, [kk-1; kk]]; % Ring segment of this trajectory (from kk-1 to kk)
                cubesats(k).tSeg = [cubesats(k).tSeg, [t(1); t(end)]]; % Time this ring segment started and ended
                cubesats(k).cost = [cubesats(k).cost, cost]; % Cost for this ring segment

                    % Update t0 and X0
                cubesats(k).X0 = X(end,1:6)';
                cubesats(k).t0 = t(end);
        
                    % Plot segment
                titleText = sprintf("Cubesat %s trajectory segment: Ring %.0f to %.0f", cubesats(k).name, kk-1, kk);
                xLabel = "Radial [m]"; yLabel = "Along Track [m]"; zLabel = "Cross Track [m]";
                if ~debug(1)
                    figNum = kk + (k-1)*length(rings);
                else
                    figNum = kk + (k-1)*length(rings) + 1;
                end
                plotSegment(cubesats(k), rings(kk), t, X, u, figNum, titleText, xLabel, yLabel, zLabel);
                
                    % Report progress
                fprintf("\tCubesat %s passed through ring %.0f in %.3f sec!\n", cubesats(k).name, kk, t(end) - t(1));
            end
            fprintf("\n\tCubesat %s finished the course in %.3f sec!\n", cubesats(k).name, cubesats(k).t(end)-cubesats(k).t(1));

                % Add trajectory to race plot
            figure(420)
            plot3(cubesats(k).X(:,1), cubesats(k).X(:,2), cubesats(k).X(:,3), '-', 'Color', cubesats(k).color, 'DisplayName', sprintf("Cubesat %s trajectory", cubesats(k).name));
        
            figure; tl = tiledlayout(1,2);
            title(tl, sprintf("Cubesat %s Adjoints", cubesats(k).name))
            nexttile
                title(sprintf("Cubesat %s Position Adjoints", cubesats(k).name))
                hold on; grid on; axis equal;
                plot3(cubesats(k).X(:,7), cubesats(k).X(:,8), cubesats(k).X(:,9), 'k.');
                view([30 35])
            nexttile
                title(sprintf("Cubesat %s Velocity Adjoints", cubesats(k).name))
                hold on; grid on; axis equal;
                plot3(cubesats(k).X(:,10), cubesats(k).X(:,11), cubesats(k).X(:,12), 'k.');
                view([30 35])
        end
        
            % Save run!
        time = datetime('now');
        runName = sprintf("Run_%s", time);
        runName = replace(runName, " ", "_");
        runName = replace(runName, ":", "");
        runName = sprintf(".\\Runs\\%s.mat", runName);
        save(runName, "cubesats", "rings", "startRing", "endRing", "courseParams", '-mat');

    case 2
        %% Choose run to analyze
        fprintf("\n\tChoose a run to analyze!\n");
        cd .\Runs\
        [file, path, ~] = uigetfile('.mat', "Pick a MAT file to analyze");
        cd ..\

        load([path,file])

        fileParts = split(file, "_");
        scenario = extractBefore(fileParts{end}, ".mat");

        menuMsg = sprintf("\n\tWould you like to save new figures?\n" + ...
                          "\t1. Yes\n" + ...
                          "\t2. No\n" + ...
                          "\tChoice (enter number of option to do): ");
        choice = input(menuMsg);

        if choice == 1
            makeFigures = true;
        else
            makeFigures = false;
        end

        %% Plot race course and example ring
            % Race course
        courseNum = 420;
        titleText = sprintf("Generated Race Course"); atEnd = false;
        xLabel = sprintf("Radial [m]"); yLabel = sprintf("Along Track [m]"); zLabel = sprintf("Cross Track [m]");
        courseFig = plotCourse(startRing, rings, endRing, cubesats, courseNum, titleText, xLabel, yLabel, zLabel, atEnd);

        if makeFigures
            angles = [0, 90, 180, 270];
            for k = angles
                view([-30 + k, 35]);
                fileName = sprintf(".\\Figures\\%s\\Course\\InitCourse_%.0fdegRot.png", scenario, k);
                saveas(courseFig, fileName);
            end
        end
        courseFig.WindowState = "minimized";

            % Example ring
        ringNum = 421;
        titleText = sprintf("Example Race Course Ring");
        xLabel = sprintf("Radial [m]"); yLabel = sprintf("Along Track [m]"); zLabel = sprintf("Cross Track [m]");
        ringFig = plotExampleRing(rings, ringNum, titleText, xLabel, yLabel, zLabel);

        if makeFigures
            for k = angles
                view([30 + k, 35]);
                fileName = sprintf(".\\Figures\\%s\\Ring\\ExampleRing_%.0fdegRot.png", scenario, k);
                saveas(courseFig, fileName);
            end
        end
        ringFig.WindowState = "minimized";

        %% Plot trajectory segments
        menuMsg = sprintf("\n\tWould you like detailed segment plots or just an overview?\n" + ...
                          "\t1. Detailed plots (one plot per segment per cubesat, including time history plots. Generates LOTS of figures)\n" + ...
                          "\t2. Overview (one plot per segment with all cubesat trajectories, no time history plots)\n" + ...
                          "\tChoice (enter number of option to do): ");
        choice = input(menuMsg);

        if choice == 1
            for k = 1:length(cubesats)
                fprintf("\nPlotting Cubesat %s's Trajectory!\n", cubesats(k).name)
                for kk = 1:length(rings)
                        % Find indices
                    kStart = find(cubesats(k).t == cubesats(k).tSeg(1,kk), 1, 'last');
                    kEnd = find(cubesats(k).t == cubesats(k).tSeg(2,kk), 1, 'first');
    
                        % Extract proper segment
                    t = cubesats(k).t(kStart:kEnd);
                    X = cubesats(k).X(kStart:kEnd, :);
                    u = cubesats(k).u(kStart:kEnd, :);
    
                        % Plot segment
                    figNum = kk + (k-1)*length(rings);
                    if kk == length(rings)
                        titleText = sprintf("Cubesat %s trajectory segment: Ring %.0f to end", cubesats(k).name, kk-1);
                    else
                        titleText = sprintf("Cubesat %s trajectory segment: Ring %.0f to %.0f", cubesats(k).name, kk-1, kk);
                    end
                    xLabel = "Radial [m]"; yLabel = "Along Track [m]"; zLabel = "Cross Track [m]";
                    figDetailed = plotSegment(cubesats(k), rings(kk), t, X, u, figNum, titleText, xLabel, yLabel, zLabel);
                        
                    if makeFigures
                        for idx = angles
                            nexttile(1);
                            view([-30 + idx, 35]);
                            fileName = sprintf(".\\Figures\\%s\\SegmentsDetailed\\DetailedSegment_%s_Ring%.0fTo%.0f_%.0fdegRot.png", scenario, cubesats(k).name, kk-1, kk, idx);
                            saveas(figDetailed, fileName);
                        end
                    end
                    figDetailed.WindowState = "minimized";

                        % Report progress
                    fprintf("\tCubesat %s passed through ring %.0f in %.3f sec!\n", cubesats(k).name, kk, cubesats(k).t(kEnd) - cubesats(k).t(kStart));
                end
                fprintf("\n\tCubesat %s finished the course in %.3f sec!\n", cubesats(k).name, cubesats(k).t(end)-cubesats(k).t(1));
    
                
            end
        else
            for k = 1:length(rings)
                    % Plot segment
                figNum = k;
                if k == length(rings)
                    titleText = sprintf("Cubesat trajectory segment: Ring %.0f to end", k-1);
                else
                    titleText = sprintf("Cubesat trajectory segment: Ring %.0f to %.0f", k-1, k);
                end
                xLabel = "Radial [m]"; yLabel = "Along Track [m]"; zLabel = "Cross Track [m]";
                figOverview = plotSegment_all(cubesats, rings(k), k, figNum, titleText, xLabel, yLabel, zLabel);

                if makeFigures
                    for idx = angles
                        nexttile(1);
                        view([-30 + idx, 35]);
                        fileName = sprintf(".\\Figures\\%s\\SegmentsOverview\\OverviewSegment_Ring%.0fTo%.0f_%.0fdegRot.png", scenario, k-1, k, idx);
                        saveas(figOverview, fileName);
                    end
                end
                figOverview.WindowState = "minimized";
            end
        end

            % Plot completed course
        courseFig = figure(courseNum);
        courseFig.WindowState = "maximized";
        view(-30, 35)
        for k = 1:length(cubesats)
            plot3(cubesats(k).X(:,1), cubesats(k).X(:,2), cubesats(k).X(:,3), '-', 'Color', cubesats(k).color, 'DisplayName', sprintf("Cubesat %s trajectory", cubesats(k).name));
        end

        if makeFigures
            for idx = angles
                view([-30 + idx, 35]);
                fileName = sprintf(".\\Figures\\%s\\Course\\FinishedCourse_%.0fdegRot.png", scenario, idx);
                saveas(courseFig, fileName);
            end
        end
        courseFig.WindowState = "minimized";
        
        %% Plot full trajectory
        for k = 1:length(cubesats)
            trajNum = 421 + k;
            titleText = sprintf("Full Trajectory of Cubesat %s", cubesats(k).name);
            trajFig = plotFullTrajectory(cubesats(k), rings, trajNum, titleText);
            if makeFigures
                fileName = sprintf(".\\Figures\\%s\\Trajectories\\Trajectory_%s.png", scenario, cubesats(k).name);
                saveas(trajFig, fileName);
            end
            trajFig.WindowState = "minimized";
        end
       

        %% Analyze race course results
        fprintf("\n\t---Analyzing Race Outcome---\n")

            % Compile costs and course times
        costs = [];
        minCosts = [];
        times = [];
        for k = 1:length(cubesats)
            cost = sum(cubesats(k).cost(1,:));
            minCost = sum(cubesats(k).cost(2,:));
            time = cubesats(k).t(end) - cubesats(k).t(1);
            costs = [costs; cost];
            minCosts = [minCosts; minCost];
            times = [times; time];
        end
        accuracies = costs - times; % cost = time cost + accuracy cost

            % Find winner (smallest cost)
        [minCost, winner] = min(costs);

            % Report results
        for k = 1:length(cubesats)
            fprintf("\nCubesat %s results: overall cost = %.3f, time = %.3f, accuracy = %.3f. \n\tMinimum possible cost: %.3f\n", cubesats(k).name, costs(k), times(k), accuracies(k), minCosts(k));
        end

        fprintf("\n\n\t~~~Cubesat %s won the race!~~~\n\n", cubesats(winner).name);
    
    case 3
        %% Choose run to animate
        fprintf("\n\tChoose a run to animate!\n");
        cd .\Runs\
        [file, path, ~] = uigetfile('.mat', "Pick a MAT file to analyze");
        cd ..\

        load([path,file])

        fileParts = split(file, "_");
        scenario = extractBefore(fileParts{end}, ".mat");

        %% Animate!
        saveMovie = true;
        movieTitle = sprintf(".\\Videos\\ASEN6020_%s_animated", scenario);
        figAnim = animateRun(cubesats, rings, saveMovie, movieTitle);

    case 4
        %% Say goodbye
        fprintf("\n\tHave a great day!\n")
        return;
end


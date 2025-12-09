%% ASEN 5254 Final Project Solution Plotter
% By: Ian Faber, 12/03/2025

%% Housekeeping
clc; clear; close all;

%% Setup
% Choose dataset
load("Run_05-May-2025_215651_Scenario2.mat");

% Read in solutions
for k = 1:4
    fprintf("Looking at Cubesat %s\n", cubesats(k).name);

    if k == 1
        solution = readmatrix("Eeny.txt");
        tSeg_1 = cumsum([0, 165.20, 88.95, 172.75, 206.50, 156.65, 129.30, 89.90, 136.50, 104.00, 128.85, 123.05, 78.35, 81.50, 34.45, 170.10, 143.05, 151.85, 87.80, 160.50, 188.05, 79.35, 153.20, 104.80, 71.90, 49.20]);
        tSeg_2 = [tSeg_1(2:end), tSeg_1(end) + 51.70];
        tSeg = [tSeg_1; tSeg_2];
    elseif k == 2
        solution = readmatrix("Meeny.txt");
        tSeg_1 = cumsum([0, 78.25, 124.50, 76.25, 18.55, 117.45, 116.45, 177.50, 180.60, 94.60, 34.45, 66.70, 111.10, 128.55, 74.85, 62.40, 85.80, 114.65, 83.05, 178.85, 181.50, 38.95, 85.00, 108.25, 156.90, 202.25]);
        tSeg_2 = [tSeg_1(2:end), tSeg_1(end) + 183.80];
        tSeg = [tSeg_1; tSeg_2];
    elseif k == 3
        solution = readmatrix("Miny.txt");
        tSeg_1 = cumsum([0, 139.50, 209.55, 113.80, 77.95, 162.40, 136.75, 44.65, 46.90, 38.20, 58.15, 43.45, 113.00, 83.45, 193.90, 105.85, 66.15, 124.85, 115.75, 66.35, 191.00, 60.80, 181.15, 83.90, 98.25, 102.40]);
        tSeg_2 = [tSeg_1(2:end), tSeg_1(end) + 88.35];
        tSeg = [tSeg_1; tSeg_2];
    else
        solution = readmatrix("Mo.txt");
        tSeg_1 = cumsum([0, 122.40, 150.35, 94.30, 81.55, 178.80, 92.05, 73.50, 182.10, 113.75, 124.50, 169.45, 154.70, 190.15, 129.70, 188.50, 88.75, 113.70, 105.85, 108.30, 150.35, 70.05, 52.95, 168.50, 85.05, 141.40]);
        tSeg_2 = [tSeg_1(2:end), tSeg_1(end) + 45.40];
        tSeg = [tSeg_1; tSeg_2];
    end

    cubesats(k).X_PVT = cubesats(k).X;
    cubesats(k).t_PVT = cubesats(k).t;
    cubesats(k).tSeg_PVT = cubesats(k).tSeg;
    cubesats(k).u_PVT = cubesats(k).u;
    cubesats(k).cost_PVT = cubesats(k).cost;
    cubesats(k).X_SST = solution(:,1:6);
    cubesats(k).t_SST = cumsum(solution(:,end));
    cubesats(k).tSeg_SST = tSeg;
    cubesats(k).u_SST = solution(:,7:9);
    cubesats(k).cost_SST = [];
    for kk = 1:length(rings)
        idx = find(cubesats(k).t_SST >= cubesats(k).tSeg_SST(2,kk), 1, 'first');
        dist = norm(cubesats(k).X_SST(idx, 1:3)' - rings(kk).center);
        time = cubesats(k).tSeg_SST(2,kk) - cubesats(k).tSeg_SST(1,kk);
        cubesats(k).cost_SST = [cubesats(k).cost_SST, [dist + time; cubesats(k).cost_PVT(2,kk)]];
    end

    fprintf("\tTotal cost for Cubesat %s: %.3f, Time taken: %.3f sec, Accuracy: %.3f m\n", cubesats(k).name, sum(cubesats(k).cost_SST(1,:)), cubesats(k).t_SST(end), sum(cubesats(k).cost_SST(1,:)) - cubesats(k).t_SST(end));
    fprintf("\tMinimum possible cost: %.3f\n", sum(cubesats(k).cost_SST(2,:)));

    cubesats(k).X = cubesats(k).X_SST;
    cubesats(k).t = cubesats(k).t_SST;
    cubesats(k).tSeg = cubesats(k).tSeg_SST;
    cubesats(k).u = cubesats(k).u_SST;
end

%% Plot race course
fig = plotCourse(startRing, rings, endRing, cubesats, 1, "Race Course", "Radial [m]", "Along-Track [m]", "Cross-Track [m]", true);

%% Plot trajectory differences
for k = 1:4
    plotFullTrajectory(cubesats(k), rings, k+1, sprintf("Cubesat %s Trajectory", cubesats(k).name));
end

return;

%% Animate race
animateRun(cubesats, rings, true, "RaceCourse_SST.mp4");





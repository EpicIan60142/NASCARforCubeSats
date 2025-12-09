%% ASEN 5254 Race Course to YAML script
% By: Ian Faber, 12/02/2025

%% Housekeeping
clc; clear; close all;

addpath("../../ASEN6020/Final Project/OfficialRuns/");

%% Setup
% Choose dataset
load("Run_05-May-2025_215651_Scenario2.mat");

% Make YAML file object
file = fopen("RaceCourse.yml",'w');

%% Populate YAML file
    % Course header
fprintf(file, "Course:\n");

    % World dimensions
fprintf(file, "  Dimensions:\n");

fig = plotCourse(startRing,rings,endRing,cubesats,1,"Race Course","Radial [m]","Along-Track [m]","Cross-Track [m]",true);

lims = zeros(3,2);
lims(1,:) = 1.1*fig.Children(3).XLim;
lims(2,:) = 1.1*fig.Children(3).YLim;
lims(3,:) = 1.1*fig.Children(3).ZLim;

for k = 1:3
    for kk = 1:2
        fprintf(file,"    - %.3f\n",lims(k,kk));
    end
end

    % Course origin dynamics
T = 180*60; % 3 hour orbit
n = 2*pi/T;
fprintf(file, "  MeanMotion: %.3e\n", n);

    % Rings
rings = [startRing; rings; endRing];
fprintf(file, "  Rings:\n");
for k = 1:length(rings)
    fprintf(file, "    ring%.0f:\n",k-1);
        % Ring center
    for kk = 1:3
        fprintf(file, "      - %.3f\n", rings(k).center(kk));
    end
        % Ring normal vector
    for kk = 1:3
        fprintf(file, "      - %.3f\n", rings(k).normal(kk));
    end
        % Ring semi-major and semi-minor axes
    for kk = 1:2
        fprintf(file, "      - %.3f\n", rings(k).S(kk,kk));
    end
        % DCM to Inertial from Ring frame
    for kk = 1:3
        for idx = 1:3
            fprintf(file, "      - %.3f\n", rings(k).NR(kk, idx));
        end
    end
end

    % Cubesats header
fprintf(file, "Cubesats:\n");

    % Cubesats
for k = 1:length(cubesats)
    fprintf(file, "  sat%.0f:\n", k-1);
        % Cubesat dynamics
    fprintf(file, "    Model: CubesatCWH\n");
        % Cubesat name, marker, and color
    fprintf(file, "    Name: %s\n", cubesats(k).name);
    fprintf(file, "    Marker: %s\n", cubesats(k).marker);
    col = split(cubesats(k).color, '#');
    fprintf(file, "    Color: %s\n", col(2));
        % Start state
    fprintf(file, "    Start:\n");
    for kk = 1:6
        fprintf(file, "      - %.3f\n", cubesats(k).X(1,kk));
    end
        % Max control
    fprintf(file, "    uMax: %.3f\n", cubesats(k).uMax);
end

fclose(file);


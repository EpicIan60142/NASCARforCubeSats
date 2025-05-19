function fig = plotFullTrajectory(cubesat, rings, figNum, titleText)
% Function that plots a cubesat's entire race trajectory
%   Inputs:
%       - cubesat: Cubesat structure for the cubesat of interest
%       - rings: Vector of ring structures for the intermediate rings of
%                the course
%       - figNum: Figure number for the trajectory figure
%       - titleText: String specifying the title for the figure
%
%   By: Ian Faber, 05/05/2025
%
    %% Make figure
fig = figure(figNum);
fig.WindowState = "maximized";
tl = tiledlayout(5,3);
title(tl, titleText);

    %% Plotting utilities
markerSize = 10;
color = cubesat.color;
ax = [];

    %% Pull out cubesat time, state, and control
t = cubesat.t;
X = cubesat.X;
u = cubesat.u;

    %% Pull out ring time, state, and control
tRings = cubesat.tSeg(2,:);
XRings = [];
normalRings = [];
velIdx = [];
for k = 1:length(rings)
    XRings = [XRings; rings(k).center'];
    normalRings = [normalRings; rings(k).normal'];
    velIdx = [velIdx; find(cubesat.t == tRings(k), 1, 'first')];
end

    %% Position
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,1), 'Color', color);
    posPlot = plot(tRings, XRings(:,1), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]"); ylabel("X [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,2), 'Color', color);
    plot(tRings, XRings(:,2), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]"); ylabel("Y [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,3), 'Color', color);
    plot(tRings, XRings(:,3), 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]"); ylabel("Z [m]");

    %% Velocity
vf = vecnorm(X(velIdx,4:6), 2, 2);
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,4), 'Color', color);
    velPlot = plot(tRings, normalRings(:,1).*vf, 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]"); ylabel("Xdot [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,5), 'Color', color);
    plot(tRings, normalRings(:,2).*vf, 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]"); ylabel("Ydot [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,6), 'Color', color);
    plot(tRings, normalRings(:,3).*vf, 'r.', 'MarkerSize', markerSize);
    xlabel("Time [sec]"); ylabel("Zdot [m/s]");

    %% Position Adjoints
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,7), 'Color', color);
    xlabel("Time [sec]"); ylabel("p_x [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,8), 'Color', color);
    xlabel("Time [sec]"); ylabel("p_y [m]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,9), 'Color', color);
    xlabel("Time [sec]"); ylabel("p_z [m]");

    %% Velocity Adjoints
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,10), 'Color', color);
    xlabel("Time [sec]"); ylabel("p_{vx} [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,11), 'Color', color);
    xlabel("Time [sec]"); ylabel("p_{vy} [m/s]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, X(:,12), 'Color', color);
    xlabel("Time [sec]"); ylabel("p_{vz} [m/s]");

    %% Control
uLim = vecnorm(u,2,2);
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,1), 'Color', color);
    plot(t, uLim, 'k--');
    plot(t, -uLim, 'k--');
    xlabel("Time [sec]");
    ylabel("u_x [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,2), 'Color', color);
    plot(t, uLim, 'k--');
    plot(t, -uLim, 'k--');
    xlabel("Time [sec]");
    ylabel("u_y [m/s^2]");
nt = nexttile; ax = [ax; nt];
    hold on; grid on;
    plot(t, u(:,3), 'Color', color);
    plot(t, uLim, 'k--');
    plot(t, -uLim, 'k--');
    xlabel("Time [sec]");
    ylabel("u_z [m/s^2]");

    %% Make legend to mark rings
lgnd = legend(posPlot, "Ring Position/Normal Vector", 'location', 'layout');
lgnd.Layout.Tile = 'east';

    %% Link axes and update figure
linkaxes(ax,'x');
drawnow;

end
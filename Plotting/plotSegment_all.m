function fig = plotSegment_all(cubesats, ring, ringIdx, figNum, titleText, xLabel, yLabel, zLabel)
% Function that plots all cubesats performance for a single segment of the 
% cubesat race course
%   Inputs:
%       - cubesats: Vector of cubesat structures for the segment of interest
%       - ring: The last ring of this segment
%       - ringIdx: Index of the current ring
%       - figNum: Figure number for this segment
%       - titleText: String specifying the title for the figure
%       - xLabel: String specifying the x axis for the 3D trajectory plot
%       - yLabel: String specifying the y axis for the 3D trajectory plot
%       - zLabel: String specifying the z axis for the 3D trajectory plot
%   Outputs:
%       - fig: Figure handle for this segment
%
%   By: Ian Faber, 05/05/2025
%

%% Create figure
fig = figure(figNum); 
fig.WindowState = "maximized";
tl = tiledlayout(1,2);
title(tl, titleText);

%% Utilities
ctrlScale = 0.25;
ax = [];

%% Plot 3D trajectory
nexttile(1)
    hold on; grid on; axis equal
        % Plot start ring and end ring
    startRing = plotRing(ring.params.lastRing, 'g-'); startRing.LineWidth = 2;
    endRing = plotRing(ring, 'r-'); endRing.LineWidth = 2;

        % Plot for each cubesat
    startAx = []; startLabels = [];
    trajAx = []; trajLabels = [];
    ctrlAx = []; ctrlLabels = [];
    for k = 1:length(cubesats)
            % Extract correct t, X, and u vectors
        kStart = find(cubesats(k).t == cubesats(k).tSeg(1,ringIdx), 1, 'last');
        kEnd = find(cubesats(k).t == cubesats(k).tSeg(2,ringIdx), 1, 'first');

            % Extract proper segment
        X = cubesats(k).X(kStart:kEnd, :);
        u = cubesats(k).u(kStart:kEnd, :);

            % Update control vector plotting indices
        uIdx = 1:40:size(u,1);

            % Update plotting utilities
        marker = cubesats(k).marker;
        color = cubesats(k).color;

            % Start position
        startAx = [startAx, scatter3(X(1,1), X(1,2), X(1,3), 'k', 'filled', 'Marker', marker)];
        
            % Plot trajectory and control vectors
        trajAx = [trajAx, plot3(X(:,1), X(:,2), X(:,3), '-', 'Color', color, 'LineWidth', 2)];
        ctrlAx = [ctrlAx, quiver3(X(uIdx,1), X(uIdx,2), X(uIdx,3), u(uIdx,1), u(uIdx,2), u(uIdx,3), ctrlScale, "filled", 'LineStyle', '-', 'Color', color)];
    
            % Update labels
        startLabels = [startLabels, sprintf("Cubesat %s Start Point", cubesats(k).name)];
        trajLabels = [trajLabels, sprintf("Cubesat %s Trajectory", cubesats(k).name)];
        ctrlLabels = [ctrlLabels, sprintf("Cubesat %s Control Vector", cubesats(k).name)];
    end

        % Labels and viewing angle
    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
    view([30 35]);

%% Make legend and update figure
    % Order legend entries properly
plotOrder = [];
plotLabels = [];
for k = 1:length(startAx)
    plotOrder = [plotOrder, startAx(k), trajAx(k), ctrlAx(k)];
    plotLabels = [plotLabels, startLabels(k), trajLabels(k), ctrlLabels(k)];
end

lgnd = legend([startRing, endRing, plotOrder], ["Start Ring", "End Ring", plotLabels], 'Location', 'layout');
lgnd.Layout.Tile = 2;

drawnow;

end
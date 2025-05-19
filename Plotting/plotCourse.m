function fig = plotCourse(startRing, rings, endRing, cubesats, figNum, titleText, xLabel, yLabel, zLabel, atEnd)
% Function that plots a generated cubesat race course
%   Inputs:
%       - startRing: Ring structure for the cubesat starting position ring
%                    of the course
%       - rings: Vector of ring structures for the intermediate rings of
%                the course
%       - endRing: Ring structure for the final cubesat rest ring
%       - cubesats: Vector of cubesat structures generated for the course
%       - figNum: Figure number for the course figure
%       - titleText: String specifying the title for the figure
%       - xLabel: String specifying the x axis for the 3D plot
%       - yLabel: String specifying the y axis for the 3D plot
%       - zLabel: String specifying the z axis for the 3D plot
%       - atEnd: Boolean indicating whether to plot Cubesats at the start
%                or end of the course
%   Outputs:
%       - fig: Course figure handle
%
%   By: Ian Faber, 05/05/2025
%

if ~exist("atEnd", 'var')
    atEnd = false;
end

    %% Make figure
fig = figure(figNum);
fig.WindowState = "maximized";
% tl = tiledlayout(1,2);
% title(tl, titleText);

    %% Plot course
% nexttile(1)
    title(titleText);
    hold on; grid on; axis equal
    
        % Plot intermediate rings
    for k = 1:length(rings)-1
        ring = scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3), 20, k, 'filled');
        normal = quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), 20, 'filled', 'k-');
        plotRing(rings(k), 'k-');
    end
        
        % Plot start and end rings
    cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
    quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), 10, 'filled', 'k-')
    plotRing(endRing, 'r-');
    
        % Plot course origin
    courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
    
    if atEnd
            % Plot cubesat ending positions
        cubeAx = []; cubeLabels = [];
        for k = 1:length(cubesats)
            cubeAx = [cubeAx, plot3(cubesats(k).X(end,1), cubesats(k).X(end,2), cubesats(k).X(end,3), 'Color', cubesats(k).color, 'Marker', cubesats(k).marker, 'MarkerFaceColor', 'k', 'MarkerSize', 5)];
            cubeLabels = [cubeLabels, sprintf("CubeSat %s", cubesats(k).name)];
        end
    else
            % Plot cubesat starting positions
        cubeAx = []; cubeLabels = [];
        for k = 1:length(cubesats)
            cubeAx = [cubeAx, plot3(cubesats(k).X(1,1), cubesats(k).X(1,2), cubesats(k).X(1,3), 'Color', cubesats(k).color, 'Marker', cubesats(k).marker, 'MarkerFaceColor', 'k', 'MarkerSize', 5)];
            cubeLabels = [cubeLabels, sprintf("CubeSat %s", cubesats(k).name)];
        end
    end
        % Labels, colorbar, and view angle
    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel); cBar = colorbar;
    cBar.Label.String = "Ring Number"; cBar.Location = "west";
    colormap("cool"); view([30 35])

    %% Legend
lgnd = legend([cubeStart, courseCenter, ring, normal, cubeAx], ["CubeSat 3\sigma Starting Sample Space", "Race Course Origin", "Course Ring", "Ring Normal Vector", cubeLabels], 'location', 'eastoutside');
% lgnd.Layout.Tile = 2;

end
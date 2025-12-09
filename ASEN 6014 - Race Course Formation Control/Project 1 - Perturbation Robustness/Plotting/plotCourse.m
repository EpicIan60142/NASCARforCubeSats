function fig = plotCourse(ringsInertial, ringsHill, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel)
% Function that plots a generated cubesat race course
%   Inputs:
%       - ringsInertial: Vector of ring structures for the intermediate rings of
%                        the course, computed in inertial space and
%                        converted to the Hill Frame
%       - ringsHill: Vector of ring structures for the intermediate rings
%                    of the course, computed entirely in the Hill frame
%       - figNum: Figure number for the course figure
%       - titleText: String specifying the title for the figure
%       - xLabel: String specifying the x axis for the 3D plot
%       - yLabel: String specifying the y axis for the 3D plot
%       - zLabel: String specifying the z axis for the 3D plot
%   Outputs:
%       - fig: Course figure handle
%
%   By: Ian Faber, 05/05/2025
%

    %% Normal vector scale
normalScale = 0.1;

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
    for k = 2:length(ringsInertial)-1
        ring = scatter3(ringsInertial(k).center(1), ringsInertial(k).center(2), ringsInertial(k).center(3), 20, k, 'filled');
        normal = quiver3(ringsInertial(k).center(1), ringsInertial(k).center(2), ringsInertial(k).center(3), ringsInertial(k).normal(1), ringsInertial(k).normal(2), ringsInertial(k).normal(3), normalScale, 'filled', 'k-');
        plotRing(ringsInertial(k), 'k-');
    end
        
        % Plot start and end rings
    startRing = ringsInertial(1);
    endRing = ringsInertial(end);
    cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
    quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), normalScale, 'filled', 'k-')
    cubeEnd = plotRing(endRing, 'r-');
    
        % Plot trajectory of each ring if it's populated
    traj1 = []; traj2 = [];
    for k = 1:length(ringsInertial)
        if ~isempty(ringsInertial(k).X)
            traj1 = plot3(ringsInertial(k).X(1,:), ringsInertial(k).X(2,:), ringsInertial(k).X(3,:), trajStyle(1));
            traj2 = plot3(ringsHill(k).X(1,:), ringsHill(k).X(2,:), ringsHill(k).X(3,:), trajStyle(2));
        end
    end
        
        % Plot course origin
    courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
    
        % Labels, colorbar, and view angle
    xlabel(xLabel); ylabel(yLabel); zlabel(zLabel); cBar = colorbar;
    cBar.Label.String = "Ring Number"; cBar.Location = "west";
    colormap("cool"); view([30 35])

    %% Legend
if ~isempty(trajStyle)
    subset = [cubeStart, cubeEnd, courseCenter, ring, normal, traj1, traj2];
    labels = ["CubeSat 3\sigma Starting Ring Initial Position", "Race Course End Ring Initial Position", "Race Course Origin", "Course Ring Initial Position", "Ring Normal Vector", trajLabel(1), trajLabel(2)];
else
    subset = [cubeStart, cubeEnd, courseCenter, ring, normal];
    labels = ["CubeSat 3\sigma Starting Ring Initial Position", "Race Course End Ring Initial Position", "Race Course Origin", "Course Ring Initial Position", "Ring Normal Vector"];
end

lgnd = legend(subset, labels, 'location', 'eastoutside');

view([-30 25]);
% lgnd.Layout.Tile = 2;

end
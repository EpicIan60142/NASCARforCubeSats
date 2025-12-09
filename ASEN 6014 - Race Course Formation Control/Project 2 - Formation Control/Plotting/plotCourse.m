function fig = plotCourse(rings, figNum, titleText, xLabel, yLabel, zLabel, trajStyle, trajLabel)
% Function that plots a generated cubesat race course
%   Inputs:
%       - rings: Vector of ring structures for the intermediate rings of
%                the course, with trajectories specified in the Hill frame
%       - figNum: Figure number for the course figure
%       - titleText: String specifying the title for the figure
%       - xLabel: String specifying the x axis for the 3D plot
%       - yLabel: String specifying the y axis for the 3D plot
%       - zLabel: String specifying the z axis for the 3D plot
%       - trajStyle: Linestyle for plotting ring trajectories
%       - trajLabel: String for labeling ring trajectories
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
    for k = 2:length(rings)-1
        ring = scatter3(rings(k).center(1), rings(k).center(2), rings(k).center(3), 20, k, 'filled');
        normal = quiver3(rings(k).center(1), rings(k).center(2), rings(k).center(3), rings(k).normal(1), rings(k).normal(2), rings(k).normal(3), normalScale, 'filled', 'k-');
        plotRing(rings(k), 'k-');
    end
        
        % Plot start and end rings
    startRing = rings(1);
    endRing = rings(end);
    cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
    quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), normalScale, 'filled', 'k-')
    cubeEnd = plotRing(endRing, 'r-');
    
        % Plot trajectory of each ring if it's populated
    traj1 = []; traj2 = [];
    for k = 1:length(rings)
        if ~isempty(rings(k).X)
            traj1 = plot3(rings(k).X(:,1), rings(k).X(:,2), rings(k).X(:,3), trajStyle(1));
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
    labels = ["CubeSat 3\sigma Starting Ring Desired Position", "Race Course End Ring Desired Position", "Race Course Origin", "Course Ring Desired Position", "Ring Normal Vector", trajLabel(1)];
else
    subset = [cubeStart, cubeEnd, courseCenter, ring, normal];
    labels = ["CubeSat 3\sigma Starting Ring Desired Position", "Race Course End Ring Desired Position", "Race Course Origin", "Course Ring Desired Position", "Ring Normal Vector"];
end

lgnd = legend(subset, labels, 'location', 'eastoutside');

view([-30 25]);
% lgnd.Layout.Tile = 2;

end
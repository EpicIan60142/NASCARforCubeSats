function fig = plotExampleRing(rings, figNum, titleText, xLabel, yLabel, zLabel)
% Function that chooses a random ring from the course and plots it
%   Inputs:
%       - rings: Vector of ring structures that make up the race course
%       - figNum: Figure number for this ring figure
%       - titleText: String specifying the title for the figure
%       - xLabel: String specifying the x axis for the 3D ring plot
%       - yLabel: String specifying the y axis for the 3D ring plot
%       - zLabel: String specifying the z axis for the 3D ring plot
%   Outputs:
%       - fig: Figure for the plotted example ring
%
%   By: Ian Faber, 05/05/2025
%
    % Choose ring
ringIdx = randi([2, length(rings)-1],1); % Choose a random intermediate ring from the list
% ringIdx = 4;
ring = rings(ringIdx);
ring.center = zeros(3,1);
normalScale = 3;
ring.S = ring.S*1000^2;

    % Make Figure
fig = figure(figNum);
fig.WindowState = 'maximized';

    % Plot ring
hold on; grid on; axis equal;
title(titleText);
ringPlot = plotRing(ring, 'r-');
ringPlot.LineWidth = 2;
normalPlot = quiver3(ring.center(1), ring.center(2), ring.center(3), ring.normal(1), ring.normal(2), ring.normal(3), normalScale, 'filled', 'b-');
xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);

ringLabel = sprintf("Course ring #%.0f - Semimajor axis = %.3f m, Semiminor axis = %.3f m", ringIdx, 1000*ring.params.a, 1000*ring.params.b);
normalLabel = sprintf("Course ring #%.0f Normal Vector - \\theta = %.3f deg, \\phi = %.3f deg", ringIdx, ring.params.theta, ring.params.phi);

legend([ringPlot, normalPlot], [ringLabel, normalLabel], 'location', 'eastoutside');

end
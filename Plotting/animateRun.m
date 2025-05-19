function fig = animateRun(cubesats, rings, saveMovie, movieTitle)
% Function that animates a race course run!
%   Inputs:
%       - cubesats: Vector of cubesat structures for the run of interest
%       - rings: Vector of ring structures for the run of interest
%       - saveMovie: Boolean indicating whether a movie of the animation is
%                    saved or not
%       - movieTitle: String specifying the title of the .MP4 to save if
%                     saveMovie is true
%   Outputs:
%       - fig: Figure handle for the animated run
%
%   By: Ian Faber, 05/06/2025
%

    %% Utilities
markerSize = 10;

    %% Find longest time vector
longestTime = 0;
longestIdx = 0;
for k = 1:length(cubesats)
    timeLength = length(cubesats(k).t);
    if timeLength > longestTime
        longestTime = timeLength;
        longestIdx = k;
    end
end

    %% Pad shorter vectors to match longest time length
for k = 1:length(cubesats)
    timeLength = length(cubesats(k).t);
    timeDiff = longestTime - timeLength;
    if timeDiff > 0
        cubesats(k).t = [cubesats(k).t; repmat(cubesats(k).t(end), timeDiff, 1)]; % Pad time vector with copies of the last time
        cubesats(k).X = [cubesats(k).X; repmat(cubesats(k).X(end,:), timeDiff, 1)];
        cubesats(k).u = [cubesats(k).u; repmat(zeros(1,3), timeDiff, 1)];
    end

end

    %% Make figure
figNum = 1;
fig = figure(figNum); fig.WindowState = "maximized";

    %% Animate the problem
movieVector = [];
dTime = 50;
for k = 1:dTime:longestTime
        % Clear figure
    clf;

        % Add tiled layout
    tl = tiledlayout(4,8);

        % Assign title
    titleText = sprintf("Race Course at t = %.3f sec", cubesats(longestIdx).t(k));
    title(tl, titleText);

        % Plot course and trajectories
    cubeAx = [];
    cubeLabels = [];
    trajAx = [];
    trajLabels = [];
    nt = nexttile(1,[4 4]);
        hold on; grid on; axis equal; 
        nt.PlotBoxAspectRatioMode = 'manual';
        nt.DataAspectRatioMode = 'manual';
        title("3D Race Course")
            % Plot intermediate rings
        for kk = 1:length(rings)-1
            ring = scatter3(rings(kk).center(1), rings(kk).center(2), rings(kk).center(3), 20, kk, 'filled');
            normal = quiver3(rings(kk).center(1), rings(kk).center(2), rings(kk).center(3), rings(kk).normal(1), rings(kk).normal(2), rings(kk).normal(3), 20, 'filled', 'k-');
            plotRing(rings(kk), 'k-');
        end
            
            % Plot start and end rings
        startRing = rings(1).params.lastRing;
        endRing = rings(end);
        cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
        quiver3(startRing.center(1), startRing.center(2), startRing.center(3), startRing.normal(1), startRing.normal(2), startRing.normal(3), 10, 'filled', 'k-')
        plotRing(endRing, 'r-');
        
            % Plot course origin
        courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
    
            % Plot current cubesat positions
        cubeAx = []; cubeLabels = [];
        for kk = 1:length(cubesats)
            cubeAx = [cubeAx, plot3(cubesats(kk).X(k,1), cubesats(kk).X(k,2), cubesats(kk).X(k,3), 'Color', cubesats(kk).color, 'Marker', cubesats(kk).marker, 'MarkerFaceColor', 'k', 'MarkerSize', markerSize)];
            cubeLabels = [cubeLabels, sprintf("CubeSat %s", cubesats(kk).name)];
        end
    
            % Plot trajectories
        for kk = 1:length(cubesats)
            trajAx = [trajAx, plot3(cubesats(kk).X(1:k,1), cubesats(kk).X(1:k,2), cubesats(kk).X(1:k,3), '-', 'Color', cubesats(kk).color)];
            trajLabels = [trajLabels, sprintf("Cubesat %s trajectory", cubesats(kk).name)];
        end

            % Labels, colorbar, and view angle
        xlabel("Radial [m]"); ylabel("Along Track [m]"); zlabel("Cross Track [m]"); 
        colormap("cool"); view(30+(k/dTime), 35)

        % Plot individual cubesats
    ctrlAx = [];
    ctrlLabels = [];
    for kk = 1:length(cubesats)
        nt = nexttile([2 2]);
            hold on; grid on; axis equal; 
            nt.PlotBoxAspectRatioMode = 'manual';
            nt.DataAspectRatioMode = 'manual';
            title(sprintf("Cubesat %s", cubesats(kk).name))
                % Cubesat position
            plot3(cubesats(kk).X(k,1), cubesats(kk).X(k,2), cubesats(kk).X(k,3), '.', 'Color', cubesats(kk).color, 'Marker', cubesats(kk).marker, 'MarkerFaceColor', 'k', 'MarkerSize', markerSize);
                % Cubesat control
            ctrlAx = [ctrlAx, quiver3(cubesats(kk).X(k,1), cubesats(kk).X(k,2), cubesats(kk).X(k,3), cubesats(kk).u(k,1), cubesats(kk).u(k,2), cubesats(kk).u(k,3), 13, 'filled', 'Color', cubesats(kk).color, 'LineWidth', 2)];
            ctrlLabels = [ctrlLabels, sprintf("Cubesat %s Control Vector", cubesats(kk).name)];

                % Make axis limits static
            x = cubesats(kk).X(k,1); y = cubesats(kk).X(k,2); z = cubesats(kk).X(k,3); 
            xlim([x-5, x+5]); ylim([y-5, y+5]); zlim([z-5, z+5]);

                % Labels and view angle
            xlabel("Radial [m]"); ylabel("Along Track [m]"); zlabel("Cross Track [m]");
            view(30, 35);
    end

        % Make legend
    plotOrder = []; 
    plotLabels = [];
    for kk = 1:length(cubesats)
        plotOrder = [plotOrder, cubeAx(kk), trajAx(kk), ctrlAx(kk)];
        plotLabels = [plotLabels, cubeLabels(kk), trajLabels(kk), ctrlLabels(kk)];
    end
    lgnd = legend([cubeStart, ring, courseCenter, plotOrder], ["Start Ring", "Course Ring", "Course Origin", plotLabels], 'location', 'layout');
    lgnd.Layout.Tile = 'east';

        % Update figure and save moviee frames
    drawnow;
    movieVector = [movieVector; getframe(fig)];

end

    % Save video if desired!
if saveMovie
    movie = VideoWriter(movieTitle, 'MPEG-4');
    movie.FrameRate = 15;

        % Open the VideoWriter object, write the movie, and close the file
    open(movie);
    writeVideo(movie, movieVector);
    close(movie);
end



end
function fig = animateRings(rings, saveMovie, movieTitle, figNum, perturb)
% Function that animates The trajectory of rings in a race course
%   Inputs:
%       - rings: Vector of ring structures for the run of interest
%       - saveMovie: Boolean indicating whether a movie of the animation is
%                    saved or not
%       - movieTitle: String specifying the title of the .MP4 to save if
%                     saveMovie is true
%       - figNum: Figure number to assign for the animation
%       - String describing the case being animated, e.g. "J2 perturbation"
%   Outputs:
%       - fig: Figure handle for the animated run
%
%   By: Ian Faber, 10/09/2025
%

    %% Utilities
markerSize = 10;

    %% Find longest time vector
longestTime = 0;
longestIdx = 0;
for k = 1:length(rings)
    timeLength = length(rings(k).t);
    if timeLength > longestTime
        longestTime = timeLength;
        longestIdx = k;
    end
end

    %% Make figure
fig = figure(figNum); fig.WindowState = "maximized";

    %% Animate the problem
movieVector = [];
dTime = 20;
for k = 1:dTime:longestTime
        % Clear figure
    clf;

        % Add tiled layout
    tl = tiledlayout(4,4);

        % Assign title
    titleText = sprintf("Race Course at t = %.3f sec, %s", rings(longestIdx).t(k), perturb);
    title(tl, titleText);

        % Plot course and trajectories
    nt = nexttile(1,[4 4]);
        hold on; grid on; axis equal; 
        nt.PlotBoxAspectRatioMode = 'manual';
        nt.DataAspectRatioMode = 'manual';
        title("3D Race Course")

            % Plot intermediate rings
        for kk = 1:length(rings)-1
            ring = scatter3(rings(kk).X(1,k), rings(kk).X(2,k), rings(kk).X(3,k), 20, kk, 'filled');
            normal = quiver3(rings(kk).X(1,k), rings(kk).X(2,k), rings(kk).X(3,k), rings(kk).normal(1), rings(kk).normal(2), rings(kk).normal(3), 0.1, 'filled', 'k-');
            rings(kk).center = [rings(kk).X(1,k); rings(kk).X(2,k); rings(kk).X(3,k)];
            plotRing(rings(kk), 'k-');
        end
            
            % Plot start and end rings
        startRing = rings(1);
        startRing.center = [startRing.X(1,k); startRing.X(2,k); startRing.X(3,k)];
        endRing = rings(end);
        endRing.center = [endRing.X(1,k); endRing.X(2,k); endRing.X(3,k)];
        cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
        quiver3(startRing.X(1,k), startRing.X(2,k), startRing.X(3,k), startRing.normal(1), startRing.normal(2), startRing.normal(3), 0.1, 'filled', 'k-')
        cubeEnd = plotRing(endRing, 'r-');
        
            % Plot course origin
        courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
    
            % Plot trajectories
        for kk = 1:length(rings)
            trajAx = plot3(rings(kk).X(1,1:k), rings(kk).X(2,1:k), rings(kk).X(3,1:k), 'b-');
            trajLabels = sprintf("Ring trajectory");
        end

            % Labels, colorbar, and view angle
        xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]"); 
        cBar = colorbar; cBar.Label.String = "Ring Number"; colormap("cool"); cBar.Location = 'westoutside';
        view(-30, 25)

        % Make legend
    lgnd = legend([cubeStart, ring, cubeEnd, courseCenter, trajAx], ["Start Ring", "Course Ring", "End Ring", "Course Origin", trajLabels], 'location', 'layout');
    lgnd.Layout.Tile = 'east';

        % Update figure and save moviee frames
    drawnow;
    movieVector = [movieVector; getframe(fig)];

end

    % Save video if desired!
if saveMovie
    movie = VideoWriter(movieTitle, 'MPEG-4');
    movie.FrameRate = 10;

        % Open the VideoWriter object, write the movie, and close the file
    open(movie);
    writeVideo(movie, movieVector);
    close(movie);
end



end
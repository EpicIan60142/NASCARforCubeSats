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
dTime = 10;
animVec = 1:dTime:longestTime;
for k = animVec
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
        for kk = 1:length(rings)
            if ~isempty(rings(kk).X)
                ring = scatter3(rings(kk).X(k,1), rings(kk).X(k,2), rings(kk).X(k,3), 20, kk, 'filled');
                normal = quiver3(rings(kk).X(k,1), rings(kk).X(k,2), rings(kk).X(k,3), rings(kk).normal(1), rings(kk).normal(2), rings(kk).normal(3), 0.1, 'filled', 'k-');
                rings(kk).center = [rings(kk).X(k,1); rings(kk).X(k,2); rings(kk).X(k,3)];
                plotRing(rings(kk), 'k-');
            end
        end
            
            % Plot start and end rings
        startRing = rings(1);
        if ~isempty(startRing.X)
            startRing.center = [startRing.X(k,1); startRing.X(k,2); startRing.X(k,3)];
            quiver3(startRing.X(k,1), startRing.X(k,2), startRing.X(k,3), startRing.normal(1), startRing.normal(2), startRing.normal(3), 0.1, 'filled', 'k-')
        end
        cubeStart = plotRing(startRing, 'g-'); cubeStart.LineWidth = 2;
        
        endRing = rings(end);
        if ~isempty(endRing.X)
            endRing.center = [endRing.X(k,1); endRing.X(k,2); endRing.X(k,3)];
        end
        cubeEnd = plotRing(endRing, 'r-'); cubeEnd.LineWidth = 2;
        
            % Plot course origin
        courseCenter = scatter3(0, 0, 0, 20, 'k', 'filled', 'h');
    
            % Plot trajectories
        for kk = 1:length(rings)
            if ~isempty(rings(kk).X)
                T = 180*60;
                cutoff = find(rings(kk).t > T, 1, 'first');
                if k < cutoff
                    start = 1;
                else
                    start = floor(k - cutoff);
                end
                trajAx = plot3(rings(kk).X(start:k,1), rings(kk).X(start:k,2), rings(kk).X(start:k,3), 'b-');
                trajLabels = sprintf("Ring trajectory");
            end
        end

            % Labels, colorbar, and view angle
        xlabel("Radial [km]"); ylabel("Along-Track [km]"); zlabel("Cross-Track [km]"); 
        cBar = colorbar; cBar.Label.String = "Ring Number"; colormap("cool"); cBar.Location = 'westoutside';
        view(45 + k*(360/longestTime), 25)

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
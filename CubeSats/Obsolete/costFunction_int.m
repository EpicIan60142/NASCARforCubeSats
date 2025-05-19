function K_int = costFunction_int(x, ring, cubesat, courseParams, opt, isFinal, debug)
% Implements the cost function for the CubeSat racing intermediate ring
% problem
%   Inputs:
%       - x: Vector of parameters to minimize the cost function with. Here,
%            that's [p0; lambda_t; lambda_v; tf]
%       - ring: Ring structure for the ring the CubeSat is trying to get
%               to, as defined in generateRing.m
%       - cubesat: Cubesat structure for the CubeSat in question, as
%                  defined by generateCubesat.m
%       - courseParams: Course parameters structure that includes the mean
%                       motion of the race course origin's orbit
%       - opt: ODE45 options via odeset
%   Outputs:
%       - K_int: Cost function value for the intermediate ring problem

    % Status indicator
persistent idx;
if isempty(idx)
    idx = 1;
end
sequence = ['/','/','/','-','-','-','\','\','\','|','|','|'];
if all(~debug)
    fprintf("\b%s",sequence(idx));
end
idx = idx + 1;
if idx > length(sequence)
    idx = 1;
end


    % Pull out p0 and tf
p0 = x(1:6);
tf = x(7);

    % Propagate [X0; p0] through CHW equations
X0 = [cubesat.X0; p0];
[t, X] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), [cubesat.t0, tf], X0, opt);

if debug(2)
        % Back out optimal control
    uOpt = [];
    for k = 1:length(t)
        pvx = X(k,10);
        pvy = X(k,11);
        pvz = X(k,12);
    
        uMax = 0;
        if length(cubesat.uMax) > 1 % Axial thrusting is at play!
            uComp = cubesat.uMax;
            if abs(pvx) > 1 % Cubesat has max x acceleration at time t
                uMax = uMax + [uComp(1); 0; 0];
            end
            if abs(pvy) > 1 % Cubesat has max y acceleration at time t
                uMax = uMax + [0; uComp(2); 0];
            end
            if abs(pvz) > 1 % Cubesat has max z acceleration at time t
                uMax = uMax + [0; 0; uComp(3)];
            end
        else
            uMax = cubesat.uMax;
        end
    
        pvHat = [pvx; pvy; pvz]/norm([pvx; pvy; pvz]);
        u = -norm(uMax)*pvHat';
        uOpt = [uOpt; u];
    end
end

    % Pull out final position
rf = X(end,1:3)';

    % Assign cost
if isFinal % Final problem - just minimize time to come to a stop
    K_int = tf - cubesat.t0;
else % Initial/intermediate problem - aim for rings
    dVec = rf - ring.center;
    K_int = norm(dVec) + (tf - cubesat.t0);
end

    % Plot iteration
if debug(2)
    % figure(69420); clf;
    %     hold on; grid on; axis equal
    %     title("Iterated Trajectory");
    %     plotRing(ring.params.lastRing, 'g-');
    %     plot3(X(:,1), X(:,2), X(:,3), 'k-');
    %     plotRing(ring, 'r-');
    %     xlabel("Radial [km]"); ylabel("Along Track [km]"); zlabel("Cross Track [km]");
    %     view([30 35]); drawnow;
    titleText = sprintf("Iterated Trajectory");
    xLabel = "Radial [m]"; yLabel = "Along Track [m]"; zLabel = "Cross Track [m]";
    plotSegment(cubesat, ring, t, X, uOpt, 69420, titleText, xLabel, yLabel, zLabel);
end

end
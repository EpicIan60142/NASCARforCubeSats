function [tOpt, XOpt, uOpt, xSolved, x0, cost] = solveTrajectory(cubesat, ring, courseParams, opt, isFinal, debug)
% Function that applies the optimal control law to navigate from one ring
% of the race course to another for the current ring problem
%   Inputs:
%       - cubesat: CubeSat currently running the course
%       - ring: Ring the CubeSat is attempting to navigate to
%       - courseParams: Course parameters structure that includes mean
%                       motion of the race course origin
%       - opt: ODE45 options via odeset
%       - isFinal: Boolean indicating whether this ring corresponds to the
%                  final problem or not
%       - debug: Boolean vector that enables or disables plotting and debug
%                features, like so:
%                - debug(1): Whether fmincon outputs a message at each
%                            iteration and plots the current guessed point
%                            true: output message and plot point
%                - debug(2): Whether the current course segment guess is
%                            plotted at each call of the cost function
%                            true: Plot current course segment guess
%                - debug(3): Whether the command window print sequence is
%                            disabled or not.
%                            true: Disable print sequence
%   Outputs:
%       - tOpt: Optimal time solved for this section of the course
%       - XOpt: Optimal CubeSat trajectory solved for this section of the 
%               course
%       - uOpt: Optimal control applied over this optimal trajectory
%       - xSolved: The resulting guess on p0 and t_f, organized as 
%                  [p0; t_f]
%       - x0: The initial guess on p0 and t_f, organized as [p0; t_f]
%       - cost: The actual and minimum possible cost for this course
%               segment, organized as [actualCost; minCost]
%
%   By: Ian Faber, 04/19/2025
%

    % Set up problem: Trying to solve for p0 and t_f such that K is minimized.
        % Set initial guesses

            % p0
r0 = cubesat.X0(1:3);
v0 = cubesat.X0(4:6);
vHat0 = v0/norm(v0);
dVec = ring.center - r0;
uMax = cubesat.uMax;
if norm(v0) == 0 % Initial ring - cubesat at rest
    p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center); 
          -dVec/norm(dVec)];
else
    p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center); 
          -(1/norm(v0))*dVec];
end

            % tf
        % Set lower bound on tf to force it to be positive and physically
        % possible given the maximum possible acceleration for this
        % cubesat. "minimum time" modeled with kinematics off of the 
        % straight line distance between the cubesat and next ring if 
        % everything lined up perfectly (it won't!). If it does, the 
        % minimum time is slightly overconfident (reduced) to ensure the 
        % optimal solution isn't skipped.
d = 0 - norm(dVec); % d = d0 - df, where d0 = 0
tMin = max((1/norm(uMax))*[-norm(v0) + sqrt(norm(v0)^2 - 2*norm(uMax)*d); -norm(v0) - sqrt(norm(v0)^2 - 2*norm(uMax)*d)]);

tf = cubesat.t0 + 10*tMin;

            % Combine initial guesses
x0 = [p0; tf];

        % Set bounds
lb = [-1e3*ones(6,1); cubesat.t0 + tMin];
ub = [1e3*ones(6,1); Inf];

    % Solve problem
if debug(1)
    options = optimoptions("fmincon", "Algorithm", "sqp", "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 3000, ...
                           "SubproblemAlgorithm", "cg", 'Display', 'iter-detailed', 'PlotFcn', 'optimplotx');
else
    options = optimoptions("fmincon", "Algorithm", "sqp", "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 3000, ...
                           "SubproblemAlgorithm", "cg", 'Display', 'none');
end

fprintf("\n") % Buffer newline for command window print sequence

[xSolved, cost] = fmincon(@(x)costFunction(x,ring,cubesat,courseParams,opt,isFinal,debug),x0,[],[],[],[],lb,ub,@(x)constraints(x,ring,cubesat,courseParams,opt,isFinal), options);

fprintf("\b") % Delete newline or last character of print sequence

    % Report constraints
[ineq, eq] = constraints(xSolved, ring, cubesat, courseParams, opt, isFinal);

if isFinal % Final problem constraints
    eqLabels = ["v_f at rest constraint", "H equivalence constraint", "H_f constraint", "X_0 constraint", "t_0 constraint"];

    [~, maxEqIdx] = max(abs(eq));

    fprintf("\t\tNo inequality constraints in final problem!\n")

    if any(maxEqIdx == 1:3)
        const = eqLabels(1);
        offset = 1;
    elseif maxEqIdx == 4
        const = eqLabels(2);
        offset = 4;
    elseif maxEqIdx == 5
        const = eqLabels(3);
        offset = 5;
    elseif any(maxEqIdx == 6:11)
        const = eqLabels(4);
        offset = 6;
    elseif maxEqIdx == 12
        const = eqLabels(5);
        offset = 12;
    else
        const = "Unknown constraint";
        offset = maxEqIdx;
    end
    fprintf("\t\tMax equality constraint magnitude: %.3e at position %.0f, %s\n", eq(maxEqIdx), maxEqIdx-offset, const);
else % Initial/intermediate problem constraints
    eqLabels = ["H equivalence constraint", "H_f constraint", "X_0 constraint", "t_0 constraint", "vHat_f pointing constraint"];
    ineqLabels = ["r_f ring distance constraint", "r_f ring plane constraint"];

    [~, maxIneqIdx] = max(abs(ineq));
    [~, maxEqIdx] = max(abs(eq));
    
    fprintf("\t\tMax inequality constraint magnitude: %.3e, %s\n", ineq(maxIneqIdx), ineqLabels(maxIneqIdx));
    
    if maxEqIdx == 1
        const = eqLabels(1);
        offset = 1;
    elseif maxEqIdx == 2
        const = eqLabels(2);
        offset = 2;
    elseif any(maxEqIdx == 3:8)
        const = eqLabels(3);
        offset = 3;
    elseif maxEqIdx == 9
        const = eqLabels(4);
        offset = 9;
    elseif any(maxEqIdx == 10:12)
        const = eqLabels(5);
        offset = 10;
    else
        const = "Unknown constraint";
        offset = maxEqIdx;
    end
    fprintf("\t\tMax equality constraint magnitude: %.3e at position %.0f, %s\n", eq(maxEqIdx), maxEqIdx-offset, const);

end

    % Propagate optimal trajectory
p0 = xSolved(1:6);
tf = xSolved(7);

X0 = [cubesat.X0; p0];

dt = 0.1;
tspan = cubesat.t0:dt:tf;
if tspan(end) ~= tf
    tspan = [tspan, tf];
end

[tOpt, XOpt] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), tspan, X0, opt);

    % Append absolute minimum cost to the cost vector for comparison.
    % Absolute minimum cost corresponds to the straight line distance time
    % taken between ring centers, assuming all velocities line up and the
    % Cubesat hits the center of the ring exactly.
cost = [cost; tMin];

    % Back out optimal control according to control law
uOpt = [];
for k = 1:length(tOpt)
    pvx = XOpt(k,10);
    pvy = XOpt(k,11);
    pvz = XOpt(k,12);

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
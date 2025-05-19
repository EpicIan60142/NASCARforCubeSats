function [tOpt, XOpt, uOpt, xSolved, x0, cost] = solveTrajectory_int(cubesat, ring, courseParams, opt, isFinal, debug)
% Function that applies the optimal control law to navigate from one ring
% of the race course to another for the intermediate ring problem
%   Inputs:
%       - cubesat: CubeSat currently running the course
%       - ring: Ring the CubeSat is attempting to navigate to
%       - courseParams: Course parameters structure that includes mean
%                       motion of the race course origin
%       - opt: ODE45 options via odeset
%   Outputs:
%       - t: Time for this section of the course
%       - X: CubeSat trajectory for this section of the course
%
%   By: Ian Faber, 04/19/2025
%

    % Set up problem: Trying to solve for p0, lambda_t, and lambda_v such that K is minimized.
        % Set initial guesses
            % Lambdas
% if cubesat.t0 == 0
%     lambda_t = -1;
%     lambda_v = 1*1e-6*ones(3,1);
%     lambda_n = 0*1e-3;
%     lambda_rf = 0*1e-3;
% else
%     % lambda_t = cubesat.optParams(7,end);
%     lambda_t = -1;
%     lambda_v = cubesat.optParams(8:10,end);
%     lambda_n = cubesat.optParams(11,end);
%     lambda_rf = cubesat.optParams(12,end);
% end

            % p0
r0 = cubesat.X0(1:3);
v0 = cubesat.X0(4:6);
vHat0 = v0/norm(v0);
dVec = ring.center - r0;
uMax = cubesat.uMax;
if norm(v0) == 0 % Initial ring - cubesat at rest
    % p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_rf*ring.params.lastRing.normal; 
    %       (eye(3) - ring.params.lastRing.normal*ring.params.lastRing.normal')*lambda_v];
    % p0 = [(1+lambda_rf)*((r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center)); 
    %       (eye(3) - ring.params.lastRing.normal*ring.params.lastRing.normal')*lambda_v];
    % p0 = [(1+lambda_rf)*(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_n*ring.params.lastRing.normal; 
    %       (eye(3) - ring.params.lastRing.normal*ring.params.lastRing.normal')*lambda_v];
    % p0 = [(1+lambda_rf)*(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_n*ring.params.lastRing.normal; 
    %       -ring.params.lastRing.normal];
    p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center); 
          -dVec/norm(dVec)];
    % p0 = ones(6,1);
else
    % p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_rf*ring.params.lastRing.normal; 
    %       (eye(3) - vHat0*vHat0')*lambda_v];
    % p0 = [(1+lambda_rf)*((r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center)); 
    %       (eye(3) - vHat0*vHat0')*lambda_v];
    % p0 = [(1+lambda_rf)*(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_n*ring.params.lastRing.normal; 
    %       (1/norm(v0))*(eye(3) - vHat0*vHat0')*lambda_v];
    % p0 = [(1+lambda_rf)*(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center) + lambda_n*ring.params.lastRing.normal; 
    %       -(1/norm(v0))*ring.params.lastRing.normal];
    p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center); 
          -(1/norm(v0))*dVec];
    % p0 = ones(6,1);
    % % p0 = [(r0 - ring.params.lastRing.center)/norm(r0 - ring.params.lastRing.center); 
    % %       -dVec/norm(dVec)];
end

            % tf
        % Set lower bound on tf to force it to be positive and physically
        % possible given the maximum possible acceleration for this
        % cubesat. "minimum time" modeled with kinematics off of the 
        % straight line distance between the cubesat and next ring if 
        % everything lined up perfectly (it won't!). If it does, the 
        % minimum time is slightly overconfident (reduced) to ensure the 
        % optimal solution isn't skipped.
% d = 0 - norm(r0 - ring.center); % d = d0 - df, where d0 = 0
d = 0 - norm(dVec); % d = d0 - df, where d0 = 0
tMin = max((1/norm(uMax))*[-norm(v0) + sqrt(norm(v0)^2 - 2*norm(uMax)*d); -norm(v0) - sqrt(norm(v0)^2 - 2*norm(uMax)*d)]);

tf = cubesat.t0 + 10*tMin;

            % Combine initial guesses
x0 = [p0; tf];

%         % Set lower bound on tf to force it to be positive and physically
%         % possible given the maximum possible acceleration for this
%         % cubesat. "minimum time" modeled with kinematics off of the 
%         % straight line distance between the cubesat and next ring if 
%         % everything lined up perfectly (it won't!). If it does, the 
%         % minimum time is slightly overconfident (reduced) to ensure the 
%         % optimal solution isn't skipped.
% % d = -ring.params.d; % d = d0 - df, df is larger
% d = 0 - norm(r0 - ring.center); % d = d0 - df, where d0 = 0
% tMin = max((1/norm(uMax))*[-norm(v0) + sqrt(norm(v0)^2 - 2*norm(uMax)*d); -norm(v0) - sqrt(norm(v0)^2 - 2*norm(uMax)*d)]);

        % Set bounds
lb = [-1e3*ones(6,1); cubesat.t0 + tMin];
ub = [1e3*ones(6,1); Inf];

%         % Define linear equalities
% sz = length(x0);
% Aeq = [
%         zeros(1,6), 1, zeros(1,sz-7)
%       ];
% beq = -1;

    % Define parameter scales
D = diag([1e-3, 1e-3, 1e-3, 1, 1, 1e-2, 1, 1, 1, 1, 1, 1, 1e-3]);

    % % Check feasibility
% optCheck = optimoptions("fmincon","Algorithm", "sqp", "EnableFeasibilityMode", true,...
%                         "SubproblemAlgorithm", "cg",'Display','iter-detailed','PlotFcn','optimplotx');
% xCheck = fmincon(@(x)0,x0,[],[],[],[],lb,ub,@(x)constraints_int(x,ring,cubesat,courseParams,opt),optCheck);

    % Solve problem
if debug(1)
    options = optimoptions("fmincon", "Algorithm", "sqp", "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 3000, ...
                           "SubproblemAlgorithm", "cg", 'Display', 'iter-detailed', 'PlotFcn', 'optimplotx');
else
    options = optimoptions("fmincon", "Algorithm", "sqp", "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 3000, ...
                           "SubproblemAlgorithm", "cg", 'Display', 'none');
end
% options = optimoptions("fmincon", "Algorithm", "sqp", "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 3000, ...
%                        "SubproblemAlgorithm", "cg", 'Display', 'final-detailed', 'PlotFcn', 'optimplotx','UseParallel',true);
% "ConstraintTolerance", 1e-1,
% options = optimoptions("fmincon", 'Display', 'iter-detailed', 'PlotFcn', 'optimplotx');
fprintf("\n")

[xSolved, cost] = fmincon(@(x)costFunction_int(x,ring,cubesat,courseParams,opt,isFinal,debug),x0,[],[],[],[],lb,ub,@(x)constraints_int(x,ring,cubesat,courseParams,opt,isFinal), options);

fprintf("\b")

    % Report constraints
[ineq, eq] = constraints_int(xSolved, ring, cubesat, courseParams, opt);

eqLabels = ["H_0 constraint", "H_f constraint", "X_0 constraint", "t_0 constraint", "vHat_f constraint", "vHat_f constraint", "r_f constraint"];
ineqLabels = ["r_f distance constraint", "r_f plane constraint", "vHat_f constraint"];

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
elseif any(maxEqIdx == 10:15)
    const = eqLabels(5);
    offset = 10;
elseif any(maxEqIdx == 16:18)
    const = eqLabels(6);
    offset = 16;
elseif maxEqIdx == 19
    const = eqLabels(7);
    offset = 19;
else
    const = "Unknown constraint";
    offset = maxEqIdx;
end
fprintf("\t\tMax equality constraint magnitude: %.3e at position %.0f, %s\n", eq(maxEqIdx), maxEqIdx-offset, const);

    % Propagate optimal trajectory
p0 = xSolved(1:6);
tf = xSolved(7);

X0 = [cubesat.X0; p0];

dt = 0.1;
% tspan = [cubesat.t0, cubesat.t0+dt:dt:tf, tf];
tspan = cubesat.t0:dt:tf;
if tspan(end) ~= tf
    tspan = [tspan, tf];
end

[tOpt, XOpt] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), tspan, X0, opt);

cost = [cost; tMin];

    % Back out optimal control
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
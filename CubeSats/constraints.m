function [c, ceq] = constraints(x, ring, cubesat, courseParams, opt, isFinal)
% Nonlinear constraint function that is passed into fmincon, places
% constraints on the trajectory and other aspects of the current ring
% problem
%   Inputs:
%       - x: Vector of parameters to minimize the cost function with. Here,
%            that's [p0; t_f]
%       - ring: Ring structure for the ring the CubeSat is trying to get
%               to, as defined in generateRing.m
%       - cubesat: Cubesat structure for the CubeSat in question, as
%                  defined by generateCubesat.m
%       - courseParams: Course parameters structure that includes the mean
%                       motion of the race course origin's orbit
%       - opt: ODE45 options via odeset
%       - isFinal: Boolean indicating whether this ring corresponds to the
%                  final problem or not
%       
%   Outputs:
%       - c: Nonlinear inequality vector, such that c(x) <= 0 for all
%            entries of c. If no inequalities, c = [].
%       - ceq: Nonlinear equality vector, such that ceq(x) == 0 for all
%              entries of c
%
%   By: Ian Faber, 04/19/2025
%

    % Pull out parameters
p0 = x(1:6);
tf = x(7);

    % Pull out ring parameters
n = ring.normal;
nHat = n/norm(n);

    % Find initial and final states for current p0 and t_f guess
X0 = [cubesat.X0; p0];
[t, X] = ode45(@(t,X)CHWEOM(t,X,cubesat,courseParams), [cubesat.t0 tf], X0, opt);

r0 = X(1,1:3)';
v0 = X(1,4:6)';

rf = X(end,1:3)';
vf = X(end,4:6)';

vfHat = vf/norm(vf);

pf = X(end,7:end)';

    % Pull out cubesat max acceleration at t0 and tf
uMax = cubesat.uMax;
uMax0 = zeros(3,1);
uMaxf = zeros(3,1);
if length(uMax) > 1  % Axial thrusting is at play!
        % t0
    if abs(dot(p0(4:6),[1;0;0])) > 0 % CubeSat has max x acceleration at t0
        uMax0 = uMax0 + [uMax(1); 0; 0];
    end
    if abs(dot(p0(4:6),[0;1;0])) > 0 % CubeSat has max y acceleration at t0
        uMax0 = uMax0 + [0; uMax(2); 0];
    end
    if abs(dot(p0(4:6),[0;0;1])) > 0 % CubeSat has max z acceleration at t0
        uMax0 = uMax0 + [0; 0; uMax(3)];
    end
        % tf
    if abs(dot(pf(4:6),[1;0;0])) > 0 % CubeSat has max x acceleration at tf
        uMaxf = uMaxf + [uMax(1); 0; 0];
    end
    if abs(dot(pf(4:6),[0;1;0])) > 0 % CubeSat has max y acceleration at tf
        uMaxf = uMaxf + [0; uMax(2); 0];
    end
    if abs(dot(pf(4:6),[0;0;1])) > 0 % CubeSat has max z acceleration at tf
        uMaxf = uMaxf + [0; 0; uMax(3)];
    end
else
    uMax0 = uMax;
    uMaxf = uMax;
end

    % Calculate initial and final gravity accelerations
g_0 = [2*courseParams.n*v0(2) + 3*courseParams.n^2*r0(1); -2*courseParams.n*v0(1); -courseParams.n^2*r0(3)];
g_f = [2*courseParams.n*vf(2) + 3*courseParams.n^2*rf(1); -2*courseParams.n*vf(1); -courseParams.n^2*rf(3)];

    % Calculate initial and final Hamiltonians
H0 = p0(1:3)'*v0 + p0(4:6)'*g_0 - norm(p0(4:6))*norm(uMax0);
Hf = pf(1:3)'*vf + pf(4:6)'*g_f - norm(pf(4:6))*norm(uMaxf);

    % Assign final constraints according to problem type
if isFinal % Final problem, want vFinal = 0 but also dynamically consistent! No rings to hit
        % Final velocity constraint
    ceq_v = vf; % vf = 0

        % Hamiltonian constraints
    ceq_Hf = Hf - (-1); % Hf = -1
    ceq_H = Hf - H0; % Hf = H0

        % Manifold constraints
    ceq_X0 = [r0; v0] - cubesat.X0;
    ceq_t0 = t(1) - cubesat.t0;

        % Assign final constraints for fmincon
    ceq = [ceq_v; ceq_H; ceq_Hf; ceq_X0; ceq_t0];
    c = [];
else % Initial and intermediate problems, want to ensure we pass through rings and keep things dynamically consistent
        % Ring inequality constraints
    c_rf_dist = norm(rf - ring.center) - min(ring.S, [], 'all'); % Want in-plane component of final intersection to be close to the ring center, i.e. |r_f - r_r,i| <= smallest size of S
    c_rf_plane = abs(dot((rf - ring.center),nHat)) - 0.1; %min(ring.S, [], 'all'); % Want normal component of final intersection to be close to the plane of the ring, i.e. dot((r_f = r_r,i), nHat) <= smallest size of S
    
        % Hamiltonian constraints
    ceq_Hf = Hf - (-1); % H_f = -1
    ceq_H = Hf - H0; % H_f = H_0
    
        % Manifold constraints
    ceq_X0 = [r0; v0] - cubesat.X0;
    ceq_t0 = t(1) - cubesat.t0;
    ceq_vHatf = vfHat - nHat;
   
        % Assign final constraints for fmincon
    ceq = [ceq_H; ceq_Hf; ceq_X0; ceq_t0; ceq_vHatf];
    c = [c_rf_dist; c_rf_plane];
end

end




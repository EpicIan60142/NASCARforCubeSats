function [t_save, X_save, u_save, doe_save, oe_d_save] = impulsiveFeedbackControl(X0, doe_r, tspan, pConst, opt)
% Function that implements an impulsive feedback control law for relative
% motion between a chief and deputy spacecraft.
%   Inputs:
%       - X0: Initial state vector in inertial coordinates organized as 
%             follows:
%             [X_c0; X_d0]
%       - doe_r: Desired orbit element differences between chief and deputy
%       - tspan: Time span to integrate over
%       - pConst: Planetary constants structure containing mu for the body
%                 of interest
%       - opt: odeoptions structure for ode45 configuration
%   Outputs:
%       - t: nx1 time vector for computed states
%       - X: Combined vector of chief and deputy states like so:
%            [X_c, X_d], where X_c, X_d are nx6 vectors
%       - u: Control effort vector like so:
%            [u_r, u_theta, u_h], where each u is an nx1 vector
%       - doe: Difference in desired and deputy orbit elements
%       - oe_d: Deputy orbit elements
%
%   By: Ian Faber, 11/18/2025
%

t_save = tspan(1);
X_save = X0';
u_save = zeros(1,3);
doe_save = [];
oe_d_save = [];

while t_save(end) < tspan(end)
    % Extract initial states
    X_c = X_save(end, 1:6)';
    X_d = X_save(end, 7:12)';
    
    % Extract chief elements
    cart.mu = pConst.mu;
    cart.rVec = X_c(1:3);
    cart.vVec = X_c(4:6);
    oe_c = convCart2ClassicOE(cart);
    
    % Extract deputy elements
    cart.rVec = X_d(1:3);
    cart.vVec = X_d(4:6);
    oe_d = convCart2ClassicOE(cart);
    M_d = convE2M(convf2E(oe_d.f, oe_d.e), oe_d.e);
    
    rHat = cart.rVec/norm(cart.rVec);
    hHat = cross(cart.rVec, cart.vVec)/norm(cross(cart.rVec, cart.vVec));
    thetaHat = cross(hHat, rHat);
    
    NH = [rHat, thetaHat, hHat];
    
    % Extract reference elements
    oe_r.a = oe_c.a + doe_r.da;
    oe_r.e = oe_c.e + doe_r.de;
    oe_r.i = oe_c.i + doe_r.di;
    oe_r.RAAN = oe_c.RAAN + doe_r.dRAAN;
    oe_r.argPeri = oe_c.argPeri + doe_r.dargPeri;
    M_c = convE2M(convf2E(oe_c.f, oe_c.e), oe_c.e);
    M_r = M_c + doe_r.dM0;
    
    % Extract difference between reference and deputy elements
    doe = [oe_r.a - oe_d.a, oe_r.e - oe_d.e, oe_r.i - oe_d.i, ...
           oe_r.RAAN - oe_d.RAAN, oe_r.argPeri - oe_d.argPeri, M_r - M_d];

    % Calculate delta V's
        % Plane change
    theta_plane = atan2(doe(4)*sin(oe_d.i), doe(3));
    f_plane = theta_plane - oe_d.argPeri;
    h = norm(cross(X_d(1:3), X_d(4:6)));
    r = norm(X_d(1:3));
    v = norm(X_d(4:6));
    dv_h = (h/r)*sqrt(doe(3)^2 + doe(4)^2*sin(oe_d.i)^2);
    phi = acos(1-((dv_h^2)/(2*v^2)));
    gamma = pi - ((pi-phi)/2);
    if f_plane < pi
        sign = 1;
    else
        sign = -1;
    end
    dV_h = [0; dv_h*cos(gamma); dv_h*sign*sin(gamma)];
    
        % Periapsis
    n_d = sqrt(pConst.mu/(oe_d.a^3));
    eta_d = sqrt(1-oe_d.e^2);
    dv_r_p = -(n_d*oe_d.a/4)*((((1+oe_d.e)^2)/eta_d)*(doe(5) + doe(4)*cos(oe_d.i)) + doe(6));
    dv_theta_p = (n_d*oe_d.a*eta_d/4)*((doe(1)/oe_d.a) + (doe(2)/(1+oe_d.e)));
    dV_p = [dv_r_p; dv_theta_p; 0];
    
        % Apoapsis
    dv_r_a = -(n_d*oe_d.a/4)*((((1-oe_d.e)^2)/eta_d)*(doe(5) + doe(4)*cos(oe_d.i)) + doe(6));
    dv_theta_a = (n_d*oe_d.a*eta_d/4)*((doe(1)/oe_d.a) - (doe(2)/(1-oe_d.e)));
    dV_a = [dv_r_a; dv_theta_a; 0];
    
    % Calculate time to each maneuver from current point
    M_plane = convE2M(convf2E(theta_plane - oe_d.argPeri, oe_d.e), oe_d.e);
    if M_plane < 0
        M_plane = M_plane + 2*pi;
    end
    M_p = 0;
    M_a = pi;
    
    if M_d < 0
        M_d_time = M_d + 2*pi;
    else
        M_d_time = M_d;
    end
    
    T = 2*pi/n_d;
    
    t_plane = (M_plane - M_d_time)/n_d;
    t_p = (M_p - M_d_time)/n_d;
    t_a = (M_a - M_d_time)/n_d;
    
    if t_p < 0
        t_p = t_p + T;
    end
    
    if t_a < 0
        t_a = t_a + T;
    end

    if t_plane < 0
        t_plane = t_plane + T;
    end
    
    % Integrate to each maneuver in sequence and apply it
    times = t_save(end) + [t_plane; t_p; t_a];
    [sortTimes, sortIdx] = sort(times);
    
    t = t_save(end);
    X = X_save(end,:);
    u = zeros(1,3);
    oe_d = [{oe_d}, t];
    dt = tspan(2) - tspan(1);
    for k = 1:length(sortTimes)
        % Define integration time
        tspan_k = t(end)+dt:dt:sortTimes(k)-dt;

        if isempty(tspan_k) || length(tspan_k) == 1
            tspan_k = t(end)+dt:dt:t(end)+2*dt;
        end
    
        % Define starting state
        X0_k = X(end,:)';
    
        % Integrate
        [t_k, X_k] = ode45(@(t,X)multiSatOrbitEOM(t,X,pConst,[struct(); struct()], false(2,1)), tspan_k, X0_k, opt);
    
        % Append integrated state to overall storage vectors
        t = [t; t_k];
        X = [X; X_k];
        u = [u; zeros(length(t_k),3)];
    
        % Apply maneuver
        switch sortIdx(k)
            case 1 % plane change happened
                u_k = dV_h;
            case 2 % periapsis happened
                u_k = dV_p;
            case 3 % apoapsis happened
                u_k = dV_a;
        end
    
        t = [t; t(end)];
        X = [X; [X(end,1:9), X(end, 10:12) + (NH*u_k)']];
        u = [u; u_k'];
    
        % Calculate resulting elements
        cart.rVec = X(end,1:3)';
        cart.vVec = X(end,4:6)';
        oe_c_k = convCart2ClassicOE(cart);

        oe_r_k.a = oe_c_k.a + doe_r.da;
        oe_r_k.e = oe_c_k.e + doe_r.de;
        oe_r_k.i = oe_c_k.i + doe_r.di;
        oe_r_k.RAAN = oe_c_k.RAAN + doe_r.dRAAN;
        oe_r_k.argPeri = oe_c_k.argPeri + doe_r.dargPeri;
        M_c_k = convE2M(convf2E(oe_c_k.f, oe_c_k.e), oe_c_k.e);
        M_r_k = M_c_k + doe_r.dM0;

        cart.rVec = X(end,7:9)';
        cart.vVec = X(end,10:12)';
        oe_d_k = convCart2ClassicOE(cart);
        M_d_k = convE2M(convf2E(oe_d_k.f, oe_d_k.e), oe_d_k.e);
        oe_d = [oe_d; {oe_d_k}, t(end)];
    
        doe_k = [oe_r.a - oe_d_k.a, oe_r.e - oe_d_k.e, oe_r.i - oe_d_k.i, ...
                 oe_r.RAAN - oe_d_k.RAAN, oe_r.argPeri - oe_d_k.argPeri, M_r_k - M_d_k];
        doe = [doe; doe_k];
    end

    t_save = [t_save; t];
    X_save = [X_save; X];
    u_save = [u_save; u];
    doe_save = [doe_save; doe];
    oe_d_save = [oe_d_save; oe_d];
end


end
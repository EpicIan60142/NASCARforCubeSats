function x_d_CL = convDeputyRL2CL(x_d_RL, r_c)
% Function that converts rectilinear deputy spacecraft coordinates to
% curvilinear deputy spacecraft coordinates
%   Inputs:
%       - x_d_RL: Deputy spacecraft rectilinear coordinates as follows:
%                 x_d_RL = [x; y; xDot; yDot]
%       - r_c: Chief orbit information vector with radius and radius rate
%              information as follows: r_c = [r; rDot]
%   Outputs:
%       - x_d_CL: Deputy spacecraft curvilinear coordinates as follows:
%                 x_d_CL = [dr; s; drDot; sDot]
%    
%   By: Ian Faber 09/14/2025
%

% Extract rectilinear coordinates and orbit info
x = x_d_RL(1);
y = x_d_RL(2);
xDot = x_d_RL(3);
yDot = x_d_RL(4);
r = r_c(1);
rDot = r_c(2);

% Map coordinates
dTheta = acot((x+r)/y);
dThetaDot = ((x+r)*yDot - y*(xDot+rDot))/((y^2) + (x+r)^2);

dr = (y/sin(dTheta)) - r;
drDot = (yDot*sin(dTheta) - y*cos(dTheta)*dThetaDot)/(sin(dTheta)^2) - rDot;

s = r*dTheta;
sDot = rDot*dTheta + r*dThetaDot;

% Assign output
x_d_CL = [dr; s; drDot; sDot];

end
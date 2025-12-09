function x_d_RL = convDeputyCL2RL(x_d_CL, r_c)
% Function that converts curvilinear deputy spacecraft coordinates to
% rectilinear deputy spacecraft coordinates
%   Inputs:
%       - x_d_CL: Deputy spacecraft curvilinear coordinates as follows:
%                 x_d_CL = [dr; s; drDot; sDot]
%       - r_c: Chief orbit information vector with radius and radius rate
%              information as follows: r_c = [r; rDot]
%   Outputs:
%       - x_d_RL: Deputy spacecraft rectilinear coordinates as follows:
%                 x_d_RL = [x; y; xDot; yDot]
%    
%   By: Ian Faber 09/14/2025
%

% Extract rectilinear coordinates and orbit info
dr = x_d_CL(1);
s = x_d_CL(2);
drDot = x_d_CL(3);
sDot = x_d_CL(4);
r = r_c(1);
rDot = r_c(2);

% Map coordinates
dTheta = s/r;
dThetaDot = (r*sDot - s*rDot)/(r^2);

x = (r+dr)*cos(dTheta) - r;
xDot = (rDot+drDot)*cos(dTheta) - (r+dr)*sin(dTheta)*dThetaDot - rDot;

y = (r+dr)*sin(dTheta);
yDot = (rDot+drDot)*sin(dTheta) + (r+dr)*cos(dTheta)*dThetaDot;


% Assign output
x_d_RL = [x; y; xDot; yDot];

end
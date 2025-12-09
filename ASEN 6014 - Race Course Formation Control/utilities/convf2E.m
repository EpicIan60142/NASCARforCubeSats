function E = convf2E(f,e)
% Function that converts eccentric anomaly to true anomaly
%   Inputs:
%       - f: True anomaly in radians
%       - e: Eccentricity
%   Outputs:
%       - E: Eccentric anomaly in radians
%
%   By: Ian Faber, 08/27/2025
%

E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));

end
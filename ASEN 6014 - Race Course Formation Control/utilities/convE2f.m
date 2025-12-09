function f = convE2f(E,e)
% Function that converts eccentric anomaly to true anomaly
%   Inputs:
%       - E: Eccentric anomaly in radians
%       - e: Eccentricity
%   Outputs:
%       - f: True anomaly in radians
%
%   By: Ian Faber, 08/27/2025
%

f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

end
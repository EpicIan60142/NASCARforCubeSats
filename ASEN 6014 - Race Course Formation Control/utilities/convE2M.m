function M = convE2M(E, e)
% Function that converts eccentric anomaly to mean anomaly
%   Inputs:
%       - E: Eccentric anomaly in radians
%       - e: Orbit eccentricity
%   Outputs:
%       - M: Mean anomaly in radians
%
%   By: Ian faber, 10/03/2025

    % Calculate mean anomaly
M = E - e*sin(E);

end
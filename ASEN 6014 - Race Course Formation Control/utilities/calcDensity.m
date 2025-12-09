function rho = calcDensity(pConst, r)
% Function that calculates atmospheric density at a given range in ECI
% coordinates
%   - Inputs:
%       - pConst: Planetary constants structure containing rho0, r0, and H
%                 for an exponential atmosphere model
%       - r: Range from the center of the planet in ECI coordinates
%   - Outputs:
%       - rho: Calculated atmospheric density at the specified range
%
%   By: Ian Faber, 02/11/2025
%
    % Extract atmospheric constants
rho0 = pConst.rho0;
r0 = pConst.r0;
H = pConst.H;

rho = rho0*exp(-(r-r0)/H);

end
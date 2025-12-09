function E = convM2E(M,e,verbose)
% Function that converts Mean Anomaly to Eccentric Anomaly using Newton's
% Method
%   Inputs:
%       - M: Mean Anomaly in radians
%       - e: Orbit eccentricity
%       - verbose: Boolean indicating whether to print results or not.
%                  Defaults to true.
%   Outputs:
%       - E: Eccentric Anomaly in radians
%
%   By: Ian Faber, 08/27/2025
%

% Set verbosity
if ~exist("verbose", "var")
    verbose = true;
end

% Define Arbitrarily small number for floating point operations to be 
% "close enough" to 0
epsilon = 1e-12; 

% Define the maximum number of iterations
maxIter = 999; 

% Define g functions for Newton's method
g = @(E) E - e*sin(E) - M;
gPrime = @(E) 1 - e*cos(E);

% Iterate E until we converge on a solution
iter = 0; % Initialize at 0 iterations
Ei = M; % Initialize guess at mean anomaly
while iter < maxIter
    E = Ei - (g(Ei)/gPrime(Ei));

    if abs(E - Ei) < epsilon
        iter = iter + 1;
        break; % Stop iteration, we have converged!
    else
        Ei = E;
        iter = iter + 1; % Keep going, still need to converge
    end
end

if verbose
    if iter < maxIter
        fprintf("\nConverged to E = %.6f rad after %.0f iterations!\n", E, iter)
    else
        fprintf("\nHit maximum iterations (%.0f) at E = %.6f rad.\n", iter, E)
    end
end

end


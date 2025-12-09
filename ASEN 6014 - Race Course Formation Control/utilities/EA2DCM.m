function mat = EA2DCM(angles, type)
% Function that converts Euler angles to a DCM for 3D rotations
%   Inputs: 
%       - angles: Vector of Euler angles corresponding to the axes in
%                 "type" in RADIANS
%       - type: Vector of axis rotations that the angles in "angles" will
%               be applied to
%   Outputs:
%       - mat: DCM for the (type) Euler Angle rotation through (angles)
%   
%   Example:
%       mat = EA2DCM(deg2rad([15, 34, 86]), [3,2,1]) 
% 
%             This will result in a (3-2-1) Euler angle rotation through 
%             the angles 15, 34, and 86 degrees. In this case, mat will 
%             be the DCM resulting from a (15) degree rotation about the 
%             (3) axis, then a (34) degree rotation about the (2) axis, 
%             then an (86) degree rotation about the (1) axis
%
%   Author: Ian Faber, 04/15/2024

M1 = @(theta) [
                1,          0,          0;
                0,  cos(theta), sin(theta);
                0,  -sin(theta), cos(theta)
              ];

M2 = @(theta) [
                cos(theta), 0,  -sin(theta);
                0,          1,            0;
                sin(theta), 0,  cos(theta)
              ];


M3 = @(theta) [
                cos(theta), sin(theta), 0;
                -sin(theta), cos(theta), 0;
                0,          0,          1
              ];

rotMats = {M1, M2, M3};

mat = (rotMats{type(3)}(angles(3)))*(rotMats{type(2)}(angles(2)))*(rotMats{type(1)}(angles(1)));

end


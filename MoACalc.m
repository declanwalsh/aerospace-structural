% function that calculates the second moment of areas for idealised structure
% B is effective boom area
% xB, yB are distances of booms from the neutral axis
function [ Ixx, Iyy, Ixy ] = MoACalc( B, yB, xB )

    narginchk(2, 3);

    if nargin == 2
        Ixx = sum(B.*(yB.^2));
    end
    if nargin == 3
        Ixx = sum(B.*(yB.^2));
        Iyy = sum(B.*(xB.^2));
        Ixy = sum(B.*xB.*yB);
    end
    
end


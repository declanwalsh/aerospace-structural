% function that calculates the basic shear flow in an idealised section
% uses the boom panel numbering provided and the general booms given
% assumes cut lies between the first two booms given (i.e. zero for first value)

% B = all boom areas of system
% yB = all boom y co-ordinates from NA of system
% numBooms = numbers of booms used (in panel defined order from before cut and defined the sign by order)
% Sy = shear force
% Ixx = 2MoA

% qBasic is returned using the order provided by numBooms (and matching that of numPanels)
function [ qBasic ] = basicSFBoxCalc( B, yB, numBooms, Sy, Ixx, xB, Sx, Iyy, Ixy)

    narginchk(5, 9)
    
    if (nargin == 5)
        BY = zeros(size(numBooms));

        for k = 2:length(numBooms)
            BY(k) = B(numBooms(k))*yB(numBooms(k)) + BY(k-1);
        end
    
        qBasic = (-Sy/Ixx)*BY;
    elseif (nargin == 9)
        BY = zeros(size(numBooms)); 
        BX = zeros(size(numBooms));
        
        for k = 2:length(numBooms)
            BY(k) = B(numBooms(k))*yB(numBooms(k)) + BY(k-1);
            BX(k) = B(numBooms(k))*xB(numBooms(k)) + BX(k-1);
        end
        
        qBasic = -(Sx*Ixx - Sy*Ixy)/(Ixx*Iyy - Ixy^2)*BX - (Sy*Iyy - Sx*Ixy)/(Ixx*Iyy - Ixy^2)*BY;
    else
        fprintf('Incorrect number of arguments (5 or 9) required')
    end
    
end


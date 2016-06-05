% Stiffener centroid function
% Author: Declan Walsh
% Last Modified: 31/05/2016

% ONLY works for horizontal or vertical components in stiffners

% INPUTS
% x, y = x and y locations of all corners (i.e. defining the section lines) 
% t = thickness of each section

% OUTPUTS
% xC, yC = x and y centroid of system
function [xC, yC] = centroidStiffener(x, y, t)

    A = zeros(size(x) - 1);
    xCIndiv = zeros(size(x) - 1);
    yCIndiv = zeros(size(y) - 1);

    for i = (1:length(x)-1)
        
        A(i) = abs((x(i+1)-x(i)).*t(i) + (y(i+1)-y(i)).*t(i)); % this line only works for horizontal/vertical lines
        xCIndiv(i) = (x(i+1)+x(i))/2;
        yCIndiv(i) = (y(i+1)+y(i))/2;
    end

    xC = sum(A.*xCIndiv)/sum(A);
    yC = sum(A.*yCIndiv)/sum(A);
    
end


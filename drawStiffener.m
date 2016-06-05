% Function that draws the stiffener defined by the user
% Author: Declan Walsh
% Last Modified: 31/05/2016

% TO DO
% highlight crippling failure corner

% INPUTS
% x, y = 
% t = thickness of each section (used to define line thickness)

function [ fig ] = drawStiffener( x, y, t, xC, yC )

    % define how much to expand the figure to capture all information
    X_SCALE = 1;
    Y_SCALE = 1;

    fig = figure;
    hold on
   
    % draw all the section lines
    for i = 1:(length(x)-1)
        line([x(i), x(i+1)], [y(i), y(i+1)], 'LineWidth', t(i));
    end
    
    % adjust the figure to capture all information
    xl = xlim;
    yl = ylim;
    
    xl(1) = xl(1)-X_SCALE;
    xl(2) = xl(2)+X_SCALE;
    yl(1) = yl(1)-Y_SCALE;
    yl(2) = yl(2)+Y_SCALE;
    
    xlim(xl)
    ylim(yl)
    
    % draw the centroid on if the information is passed into the function
    if nargin == 5
       plot(xC, yC, 'r*') 
    end
    
end
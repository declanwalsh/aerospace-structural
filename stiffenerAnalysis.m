% Stiffener analysis function
% Author: Declan Walsh
% Last Modified: 31/05/2016

% uses the thin walled approximation to analyse stiffeners
% all stiffeners are composed of a base/web and flange

% INPUTS
% stiffener = structure containing all stiffener geometry
% draw = true/false value to draw the stiffener and resultant info

% OUTPUTS
% x, y = x and y positions of each corner (defining each section line)
% t = thickness of each section

function [ x, y, t ] = stiffenerAnalysis( stiffener, draw )

    % analysis depends on the type of stiffener
    % same parameters used in different ways for each case
    switch stiffener.type
        case 'U'
            % defined by 6 points (4 corners and 2 ends)
            % forms 5 lines
            x = [0, stiffener.base_l, stiffener.base_l, stiffener.base_l + stiffener.flange_l, ...
                stiffener.base_l + stiffener.flange_l, stiffener.base_l + stiffener.flange_l + stiffener.base_l];
            y = [0, 0, stiffener.web_l, stiffener.web_l, 0, 0];
            t = [stiffener.t_base, stiffener.t_web, stiffener.t_flange, stiffener.t_web, stiffener.t_base];
        case 'Z'
            % defined by 4 points (2 corner and 2 ends)
            % forms 3 lines
            x = [0, stiffener.base_l, stiffener.base_l, stiffener.base_l + stiffener.flange_l];
            y = [0, 0, stiffener.web_l, stiffener.web_l];
            t = [stiffener.t_base, stiffener.t_web, stiffener.t_flange];
        case 'T'
            % defined by 5 points (treated as two back to back L)
            % forms 4 lines (middle 2 are 2 halfs of the web)
            x = [0, stiffener.base_l, stiffener.base_l, stiffener.base_l,...
                stiffener.base_l + stiffener.base_l];
            y = [0, 0, stiffener.web_l, 0, 0];
            t = [stiffener.t_base, stiffener.t_web/2, stiffener.t_web/2, stiffener.t_base];
        otherwise
            fprintf('Unknown stiffener type encountered');
    end

    % calculate centroids
    [xC, yC] = centroidStiffener(x, y, t)

    % draw the stiffeners if desired
    if (draw == 1)
        drawStiffener(x, y, t, xC, yC);
    end
    
    % calculate the second moment of areas of the stiffner
   % I = 

end


% Function that analyses wing panel parameters
% Author: Declan Walsh
% Last Modified: 31/05/2016

% INPUTS 
% skin, stiffener = structures containing various geometrical parameters

% OUTPUT 
% F = maximum stresses for various failure modes

function [F] = panel(skin, stiffener)

% CONTROLLABLE PARAMETERS
% t_skin = skin thickness (mm) - affects panel buckling, 
% pitch_stiffener = stiffener pitch (mm) - affects panel buckling, 

% t_stiffener = thickness of stiffener flange/web (mm) - affects stiffener local buckling
% b_stiffener = length of stiffener flange/web (mm) -affects stiffener local buckling

% Aluminium values
E = 70e9;
v = 0.33;

k_flat_4S = 4; % buckling coefficient for flat plate simply supported on all 4 sides
k_flat_1F3S = 0.43; % buckling coefficient for flat plate simply supported on 3 sides and free on one

F.panel_buckling = k_flat_4S*pi^2*E/(12*(1-v^2))*(skin.t/skin.pitch_stiffener)^2;
% F.stiffener_local_buckling = k_flat_1F3S*pi^2*E/(12*(1-v^2))*(stiffener.t/stiffener.b_stiffener);


end


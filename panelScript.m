% script to test the analysis functions
% Author: Declan Walsh
% Last Modified: 31/05/2016

close all
clear all

skin.t = 1;
skin.pitch_stiffener = 50;

stiffener.type = 'Z';
stiffener.base_l = 5;
stiffener.web_l = 2;
stiffener.flange_l = 5;
stiffener.t_base = 2;
stiffener.t_web = 4;
stiffener.t_flange = 2;

[x, y, t] = stiffenerAnalysis(stiffener, 1)

F = panel(skin, stiffener);
% Script to solve loading of bolts in fastener groups
% Q Joints - 2
% AERO3410 Tutorials
% all dimensions are given from a common point

clear all
close all

%% BOLT PARAMETERS

D = 0.5; % diameter of bolts
t = 0.05; % thickness of fastener skin

A = pi*D^2/4; % area of bolts (assumes circular)

Px = 0; % direct force on bolts in x direction
Py = -400; % direct force on bolts in y direction

xLoading = 0; % distance of applied load from origin in x-component
yLoading = 0; % distance of applied load from origin in y-component

Fshear = 40000; % fastener shear failure stress (psi)
Psu = Fshear*A; % critical shear fastener failure load

Fbr = 82000; % bearing failure stress (psi)
Pbru = Fbr*t*D; % critical bearing stress load

%% GEOMETRY OF FASTENER GROUPS

clip1 = [1, 2];
clip1A = [A, A];
clip1x = [0.475, 1.005];
clip1y = [0, 0.53];
numBolts1 = length(clip1);

centroid1x = sum(clip1A.*clip1x)/sum(clip1A);
centroid1y = sum(clip1A.*clip1y)/sum(clip1A);

r1 = sqrt((clip1x - centroid1x).^2 + (clip1y - centroid1y).^2);

clip2 = [3, 4];
clip2A = [A, A];
clip2x = [0.425, 1.175];
clip2y = [0, 0];
numBolts2 = length(clip2);

centroid2x = sum(clip2A.*clip2x)/sum(clip2A);
centroid2y = sum(clip2A.*clip2y)/sum(clip2A);

r2 = sqrt((clip2x - centroid2x).^2 + (clip2y - centroid2y).^2);

%% LOADING OF FASTENERS

M1 = Px*(yLoading - centroid1y) + Py*(xLoading - centroid1x);
M1 = 0;

clip1DirectLoadX = Px/numBolts1;
clip1DirectLoadY = Py/numBolts2;

clip1MomentLoadX = M1.*clip1A.^2.*(clip1y - centroid1y)/sum(clip1A.^2.*r1.^2);
clip1MomentLoadY = M1.*clip1A.^2.*(clip1x - centroid1x)/sum(clip2A.^2.*r2.^2);

clip1LoadX = clip1DirectLoadX + clip1MomentLoadX;
clip1LoadY = clip1DirectLoadY + clip1MomentLoadY;

clip1Load = sqrt(clip1LoadX.^2 + clip1LoadY.^2);

%clip2LoadX = 
%clip2LoadY = 

%% FAILURE

MoSFastenerShear = Psu/max(clip1Load) - 1;
fprintf('The margin of safety on fastener shear failure is: %.3f\n', MoSFastenerShear);

MoSBearing = Pbru/max(clip1Load) - 1;
fprintf('The margin of safety on bearing failure is: %.3f\n', MoSBearing);
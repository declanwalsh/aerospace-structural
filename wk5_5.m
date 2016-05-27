%%
% ROTATION OF AXIS/FIBRES
clear all;

theta = -7.5; % orientation of fibres relative to axis

c = cosd(theta);
s = sind(theta);

rot = [c^2, s^2, 2*c*s; s^2, c^2, -2*c*s; -c*s, c*s, c^2-s^2]; % rotation matrix
rotInv = [c^2, s^2, -2*c*s; s^2, c^2, 2*c*s; c*s, -c*s, c^2-s^2]; % inverse of rotation matrix equivalent to rot(-theta)

%%
% STRESS IN LOCAL FIBRE CO-ORDINATE SYSTEM

Sx = 500;
Sy = 20;
Sxy = 45;

sInitial = [Sx; Sy; Sxy]
sFibre = rot*sInitial
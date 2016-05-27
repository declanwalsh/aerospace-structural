%%
% MATERIAL PROPERTIES
clear all
clc

E11 = 150;
E22 = 7.5;
G12 = 3.8;
v12 = 0.25;
v21 = E22*v12/E11;

XT = 450;
XC = 350;
YT = 30;
YC = 40;
Smax = 35;

theta = (-22.5+5*3); % angle of rotation in degrees
Sx = 500; % stress in x direction in MPa
Sy = 20; % stress in y direction in MPa
Sxy = (-45 + 10*9); % shear stress in MPa

S = [Sx; Sy; Sxy]

%%
% REDUCED STIFFNESS COEFFICENT MATRIX
Q = zeros(3,3);
Q(1,1) = E11/(1-v12*v21);
Q(1,2) = v21*E11/(1-v12*v21);
Q(2,1) = v12*E22/(1-v12*v21);
Q(2,2) = E22/(1-v12*v21);
Q(3,3) = G12;

Q

%%
% ROTATION OF AXIS/FIBRES

c = cosd(theta);
s = sind(theta);

rot = [c^2, s^2, 2*c*s; s^2, c^2, -2*c*s; -c*s, c*s, c^2-s^2]; % rotation matrix
rotInv = [c^2, s^2, -2*c*s; s^2, c^2, 2*c*s; c*s, -c*s, c^2-s^2]; % inverse of rotation matrix equivalent to rot(-theta)

Qrot = rotInv*Q*rot

Sprin = rot*S
%% FAILURE THEORIES
% Maximum Stress

if(Sprin(1) > XT || Sprin(1) < XC || Sprin(2) > YT || Sprin(2) < YC || Sprin(3) > Smax)
    fprintf('FAILURE - Maximum stress failure theory criteria met\n')
else
    fprintf('Maximum stress failure theory does not predict failure\n')
end

tsaiHill = (Sprin(1)/XT)^2 + (Sprin(2)/YT)^2 - (Sprin(1)*Sprin(2)/(YT*XT)) + (Sprin(3)/Smax)^2
if tsaiHill >= 1
    fprintf('FAILURE - Tasi Hill theory criteria met\n')
else
    fprintf('Tsai Hill failure theory does not predict failure\n')
end

F1 = 1/XT + 1/XC;
F2 = 1/YT + 1/YC;
F3 = 0; % for plane stress
% F4, F5 and F6 = 0

F11 = -1/(XT*XC);
F22 = -1/(YT*YC);
F66 = 1/(Smax^2);
F12 = 0;

tsaiWu = F1*Sprin(1) + F2*Sprin(2) + F11*Sprin(1)^2 + F12*Sprin(1)*Sprin(2)  + F22*Sprin(2)^2 + F66*Sprin(3)^2
if tsaiWu >= 1
    fprintf('FAILURE - Tsai Wu theory criteria met\n')
else
    fprintf('Tsai Wu failure theory does not predict failure\n')
end

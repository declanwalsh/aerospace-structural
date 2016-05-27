%%
% MATERIAL PROPERTIES
clear all
clc

NUM_LAYERS = 5;
theta = -45; % orientation of fibres relative to axis

% reduced stiffness coefficient of ply 0/90
Qxx = 74;
Qyy = 74;
Qxy = 2.9;
Qss = 7.2;

Qnorm = [Qxx, Qxy, 0; Qxy, Qyy, 0; 0, 0, Qss]

%%
Qinv = inv(Qnorm);
E1 = 1/Qinv(1,1)
E2 = 1/Qinv(2,2)
v12 = -Qinv(2,1)*E1
v21 = -Qinv(1,2)*E2
G12 = 1/Qinv(3,3)

%%
% ROTATION OF AXIS/FIBRES

c = cosd(theta);
s = sind(theta);

rot = [c^2, s^2, 2*c*s; s^2, c^2, -2*c*s; -c*s, c*s, c^2-s^2]; % rotation matrix
rotInv = [c^2, s^2, -2*c*s; s^2, c^2, 2*c*s; c*s, -c*s, c^2-s^2]; % inverse of rotation matrix equivalent to rot(-theta)

Qrot = rotInv*Qnorm*rot
Qcore = zeros(3,3);

h = [-55, -30, -5, 5, 30, 55];
Q = cat(3, Qrot, Qnorm, Qcore, Qnorm, Qrot);

A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for i = 1:NUM_LAYERS
    A = A + Q(:, :, i).*(h(i+1) - h(i));
    B = B + Q(:, :, i).*(h(i+1)^2 - h(i)^2);
    D = D + Q(:, :, i).*(h(i+1)^3 - h(i)^3);
end

B = B/2;
D = D/3;

% ensures the entire B matrix is approimately 0
assert(all(all(B < 1e-6)), 'B should be equal to zero for a symmetric laminate');
ABBD = [A, B; B, D]
abcd = inv(ABBD)
dxx = abcd(4,4)

% bending stiffness
b = 20; % cross sectional width in mm
EI = b/dxx

P = 500; % load in Newtons
b = 200; % distance from load to support in mm
L = 500; % total beam length in mm
centralDeformation = P*b*(3*L^2 - 4*b^2)/(48*EI) % central deformation in mm
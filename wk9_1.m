% script to solve an idealised 2 cell wing in pure single direction loading of shear
%

B = [1290, 1936, 645, 645, 1936, 1290]; % effective boom area -> 1, 2, etc.
yB = [-127, -203, -101, 101, 203, 127]; % boom y locations

panel1Length = [647, 254, 647, 406]; % panels of box 1 (56, 61, 12, 25) lengths
panel1MomentArm = [0, -635, 0 , 0];
panel1Boom = [5, 6, 1, 2]; % booms associated with each panel (in same order)
panel1Thickness = [0.915, 1.625, 0.915, 2.032]; % panel 1 thickness

panel2Length = [775, 406, 775, 202]; % panels of box 2 (45, 52, 23, 34) lengths
panel2MomentArm = [0, 0, 0, 763 ];
panel2Boom = [4, 5, 2, 3];
panel2Thickness = [0.559, 2.032, 0.559, 1.22]; % panel 2 thicknesses

panel1Delta = panel1Length./panel1Thickness;
panel2Delta = panel2Length./panel2Thickness;

A1 = 232000; % area of box  1
A2 = 258000; % area of box 2

G = 1; % placeholder shear modulus

Sy = 44500; % vertical load applied
Ixx = MoACalc(B, yB); % calculate 2MoA

qB1 = basicSFBoxCalc( B, yB, panel1Boom, Sy, Ixx);
qB2 = basicSFBoxCalc( B, yB, panel2Boom, Sy, Ixx);

% find common booms and panel indexes relating to the booms
[commonB, panel1CommonIdx, panel2CommonIdx] = commonBooms(panel1Boom, panel2Boom);

% rate of twist equation (Q11*qs01 + Q12*qs02 = Q13)
Q11 = 1/(2*A1*G)*sum(panel1Delta) + 1/(2*A2*G)*panel2Delta(panel2CommonIdx);
Q12 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) - 1/(2*A2*G)*sum(panel2Delta);
Q13 = (1/(2*A2*G))*sum(qB2.*panel2Delta) - (1/(2*A1*G))*sum(qB1.*panel1Delta);

% moment equilibrium (Q21*qs01 + Q22*qs02 = Q23)
Q21 = 2*A1;
Q22 = 2*A2;
Q23 = sum(qB1.*panel1Length.*panel1MomentArm) - sum(qB2.*panel2Length.*panel2MomentArm);

A = [Q11, Q12; Q21, Q22];
B = [Q13; Q23];

qs0 = inv(A)*B;

qTotal1 = qB1 + qs0(1)
qTotal2 = qB2 + qs0(2)
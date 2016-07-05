function [ FoS_incTorsion, FoSUlt_incTorsion, FoSSkin, FoSSkin_curv, FoSir, FoScr, FoScripPure, FoScripAdjusted, FoSStringer ] = analysisTotalLegacy( LOAD_FACTOR, SAVE_DATA )

% function to solve an idealised 3 cell wing in combination of shear and moment
% Author: Declan Walsh
% For AERO3410 wing analysis assignment
% all analysis has been performed in Imperial Units

ROUGH_RELATIVE_ERROR = 0.05; % used as percentage margin to check results
PRECISE_RELATIVE_ERROR = 1e-10;

%% BASIC GEOMETRY AND CONSTANTS

PA_TO_PSI = 0.000145038; % conversion unit for Pascals to psi
NMM_TO_PI = 5.71015; % conversion unit for N/mm to lb/in
FT_TO_IN = 12; % conversion unit for feet to inches
MPH_TO_FPS = 17.6; % conversion unit for mph to feet per second

KILO = 1000;

G = 28e9 * PA_TO_PSI; % shear modulus of Al-2024 sheet metal clad (Pa)
F_YIELD = 310e6 * PA_TO_PSI; % yield stress of Al-2024 sheet metal clad (ksi)
F_ULT = 430e6 * PA_TO_PSI;

CHORD = 6*FT_TO_IN; % length of chord of Beechcraft King Air B200 (in)

t_skin = 0.1; % skin thickness (in) - controls skin buckling
t_stiffner = 0.3; % stiffner thickness (in)

STRINGER_WIDTH = 0.5; % (in)
STRINGER_DEPTH = 1; % (in) - controls stiffener column buckling
t_stringer = 0.2;

% aircraft data
% http://www.mrmoo.net/pilot/BE200/Older%20CRH%20split%20up/03%20Limitations.pdf
BANK_ANGLE_MAX = acosd(1/LOAD_FACTOR); % bank angle of aircraft that gives the load factor
MTOW = 12500; % maximum take off weight (lb)
WEIGHT_WING = 2000; % weight per wing (lb)
SPAN = 60*FT_TO_IN; % span (in)
EQ_MTOW = LOAD_FACTOR * MTOW;
EQ_WEIGHT_WING = LOAD_FACTOR*WEIGHT_WING;

% drag calc data
DENSITY = 0.0765/(FT_TO_IN^3); % density of air at 59F in lbm/inch^3
g = 32.174*FT_TO_IN; % acceleration due to gravity in in per second^2
V_MAX = 339*MPH_TO_FPS;
AREA_WING = 303*FT_TO_IN*FT_TO_IN;
q_MAX = 0.5*DENSITY/g*V_MAX^2;
CD0 = 0.02;
K = 1/(pi*9.8*0.85);
CD = CD0 + K*((EQ_MTOW/(AREA_WING*q_MAX))^2);
D = q_MAX*AREA_WING*CD;

SF_D = D;
BM_D = 0.5 * D * (SPAN/4);

% load data
SF = 0.5 * EQ_MTOW - EQ_WEIGHT_WING;
BM = 0.4*(SPAN/2)*SF;

% shear forces act at the quarter chord and through the centreline
SF_xLocation = CHORD/4;
SF_yLocation = 0;

% airfoil modelled off NACA-23018 - values fixed by eye
% lengths of stiffners (in terms of chord)
STIFFNER_1 = 0.14*CHORD;
STIFFNER_2 = 0.18*CHORD;
STIFFNER_3 = 0.1*CHORD;

% locations of stiffners from leading edge (in terms of chord)
LENGTH_1 = 0.15*CHORD;
LENGTH_2 = 0.325*CHORD;
LENGTH_3 = 0.325*CHORD;

frontEllipseApprox = pi*sqrt((LENGTH_1^2 + (STIFFNER_1/2)^2)/2); % approximation to ellipse length
frontEllipseLength = (pi/2)*(3*(LENGTH_1 + STIFFNER_1/2) - sqrt((3*LENGTH_1+STIFFNER_1/2)*(LENGTH_1 + 3* STIFFNER_1/2))); % more accurate ellipse length
assert(abs(frontEllipseApprox/frontEllipseLength - 1) < ROUGH_RELATIVE_ERROR, 'Ellipse length equations may be inaccurate - two methods have diverging results');

%% IDEALISATION OF STRUCTURE

% boom y locations
yB = [linspace(-STIFFNER_1/2, -STIFFNER_2/2, 5), ... % booms 1, 2, 3, 4, 5
    linspace(-STIFFNER_2/2, -STIFFNER_3/2, 5), ... % booms 5, 6, 7, 8, 9
    linspace(STIFFNER_3/2, STIFFNER_2/2, 5) ... % booms 10, 11, 12, 13, 14
    linspace(STIFFNER_2/2, STIFFNER_1/2, 5)]; % booms 14, 15, 16, 17, 18

xB = [linspace(LENGTH_1, LENGTH_1 + LENGTH_2, 5), ...
    linspace(LENGTH_1 + LENGTH_2, LENGTH_1 + LENGTH_2 + LENGTH_3, 5)];

yB = yB([1, diff(yB)]~=0); % remove duplicates occuring due to doubling up
xB = xB([1, diff(xB)]~=0);

xB = [xB, fliplr(xB)]; % top half is equivlant to bottom half

% panel 1 defined (elliptical - leading edge)
panel1Length = [frontEllipseLength, STIFFNER_1]; % panels of box 1 (18-1 outer, 1-18 inner) lengths
panel1Boom = [18, 1]; % booms associated with each panel (in same order)
panel1Thickness = [t_skin, t_stiffner]; % panel 1 thickness (outer, inner)
panel1Y = yB(panel1Boom);

%panel1MomentArm = mean([panel1Y(1:end-1);panel1Y(2:end)]);
%panel1MomentArm(end + 1) = (panel1Y(end) + panel1Y(1))/2;
panel1MomentArm(1) = -SF_xLocation;
panel1MomentArm(2) = -SF_xLocation + LENGTH_1;

% panel 2 defined (trapezoidal - centre)
INDIV_LENGTH_2 = sqrt(LENGTH_2^2 + ((STIFFNER_1 - STIFFNER_2)/2)^2)/4; % length of panels between longerons
panel2Length = [repmat(INDIV_LENGTH_2 ,1 , 4), STIFFNER_2, repmat(INDIV_LENGTH_2 , 1, 4), STIFFNER_1]; % panels of box 2 (14-15, ... 17-18, 18-1, 1-2, ..., 4-5) lengths
panel2Boom = [18, 17, 16, 15, 14, 5, 4, 3, 2, 1];
panel2Thickness = [repmat(t_skin, 1, 4), t_stiffner, repmat(t_skin, 1, 4), t_stiffner]; % panel 2 thicknesses
panel2Y = yB(panel2Boom);

panel2MomentArm = mean([panel2Y(1:end-1);panel2Y(2:end)]);
panel2MomentArm(5) = SF_xLocation - LENGTH_1;
panel2MomentArm(end + 1) = SF_xLocation - LENGTH_1 - LENGTH_2;

% panel 3 defined (trapezoidal - trailing edge)
INDIV_LENGTH_3 = sqrt(LENGTH_3^2 + ((STIFFNER_2 - STIFFNER_3)/2)^2)/4; % length of panels between stringers
panel3Length = [repmat(INDIV_LENGTH_3 ,1 , 4), STIFFNER_2, repmat(INDIV_LENGTH_3 , 1, 4), STIFFNER_3]; % panels of box 2 (14-5, ... 17-18, 18-1, 1-2, ..., 4-5) lengths
panel3Boom = [10, 11, 12, 13, 14, 5, 6, 7, 8, 9];
panel3Thickness = [repmat(t_skin, 1, 4), t_stiffner, repmat(t_skin, 1, 4), t_stiffner]; % panel 2 thicknesses
panel3Y = yB(panel3Boom);

panel3MomentArm = mean([panel3Y(1:end-1);panel3Y(2:end)]);
panel3MomentArm(end + 1) = (panel3Y(end) + panel3Y(1))/2;
panel3MomentArm(5) = -SF_xLocation + LENGTH_1 + LENGTH_2;
panel3MomentArm(end) = -SF_xLocation + LENGTH_1 + LENGTH_2 + LENGTH_3; 

panel1Delta = panel1Length./panel1Thickness;
panel2Delta = panel2Length./panel2Thickness;
panel3Delta = panel3Length./panel3Thickness;

A1 = pi*LENGTH_1*STIFFNER_1/2; % area of box  1 (ellipse)
A2 = ((STIFFNER_1+STIFFNER_2)/2)*LENGTH_2; % area of box 2 (trapezoid)
A3 = ((STIFFNER_2+STIFFNER_3)/2)*LENGTH_3; % area of box 3 (trapezoid);

%% BOOM CALCULATIONS

A_stringer = STRINGER_DEPTH*t_stringer + (STRINGER_WIDTH-t_stringer)*t_stringer + (STRINGER_WIDTH-t_stringer)*t_stringer;
A_spar_cap = 2*A_stringer; % spar caps assumed to be back to back stringers for analysis

A_actual = [A_spar_cap, A_stringer, A_stringer, A_stringer, A_spar_cap, A_stringer, A_stringer, A_stringer, A_stringer, ...
    A_stringer, A_stringer, A_stringer, A_stringer, A_spar_cap, A_stringer, A_stringer, A_stringer, A_spar_cap];

% assuming a pure bending moment to calculate relative stresses
% assuming symmetry and only solving for the bottom half
A_effective(1) = A_spar_cap + panel1Length(1)*panel1Thickness(1)*(2+yB(18)/yB(1))/6 ...
    + panel1Length(2)*panel1Thickness(2)*(2+yB(18)/yB(1))/6 ...
    + panel2Length(9)*panel2Thickness(9)*(2+yB(2)/yB(1))/6;

A_effective(2) = A_stringer + panel2Length(9)*panel2Thickness(9)*(2+yB(1)/yB(2))/6 + panel2Length(8)*panel2Thickness(8)*(2+yB(3)/yB(2))/6;
A_effective(3) = A_stringer + panel2Length(8)*panel2Thickness(8)*(2+yB(2)/yB(3))/6 + panel2Length(7)*panel2Thickness(7)*(2+yB(4)/yB(3))/6;
A_effective(4) = A_stringer + panel2Length(7)*panel2Thickness(7)*(2+yB(3)/yB(4))/6 + panel2Length(6)*panel2Thickness(6)*(2+yB(5)/yB(4))/6;

A_effective(5) = A_spar_cap + panel2Length(6)*panel2Thickness(6)*(2+yB(4)/yB(5))/6 ...
    + panel2Length(5)*panel2Thickness(5)*(2+yB(14)/yB(5))/6 ...
    + panel3Length(6)*panel3Thickness(6)*(2+yB(6)/yB(5))/6;

A_effective(6) = A_stringer + panel3Length(6)*panel3Thickness(6)*(2+yB(5)/yB(6))/6 + panel3Length(7)*panel3Thickness(7)*(2+yB(7)/yB(6))/6;
A_effective(7) = A_stringer + panel3Length(7)*panel3Thickness(7)*(2+yB(6)/yB(7))/6 + panel3Length(8)*panel3Thickness(8)*(2+yB(8)/yB(7))/6;
A_effective(8) = A_stringer + panel3Length(8)*panel3Thickness(8)*(2+yB(7)/yB(8))/6 + panel3Length(9)*panel3Thickness(9)*(2+yB(9)/yB(8))/6;

A_effective(9) = A_spar_cap + panel3Length(9)*panel3Thickness(9)*(2+yB(8)/yB(9))/6 ...
    + panel3Length(10)*panel3Thickness(10)*(2+yB(10)/yB(9))/6;

%% ALTERNATE BOOM CALC

% assuming a pure bending moment to calculate relative stresses
% assuming symmetry and only solving for the bottom half
% A_effectiveD(1) = A_spar_cap + panel1Length(1)*panel1Thickness(1)*(2+xB(18)/xB(1))/6 ...
%     + panel1Length(2)*panel1Thickness(2)*(2+xB(18)/xB(1))/6 ...
%     + panel2Length(9)*panel2Thickness(9)*(2+xB(2)/xB(1))/6;
% 
% A_effectiveD(2) = A_stringer + panel2Length(9)*panel2Thickness(9)*(2+xB(1)/xB(2))/6 + panel2Length(8)*panel2Thickness(8)*(2+xB(3)/xB(2))/6;
% A_effectiveD(3) = A_stringer + panel2Length(8)*panel2Thickness(8)*(2+xB(2)/xB(3))/6 + panel2Length(7)*panel2Thickness(7)*(2+xB(4)/xB(3))/6;
% A_effectiveD(4) = A_stringer + panel2Length(7)*panel2Thickness(7)*(2+xB(3)/xB(4))/6 + panel2Length(6)*panel2Thickness(6)*(2+xB(5)/xB(4))/6;
% 
% A_effectiveD(5) = A_spar_cap + panel2Length(6)*panel2Thickness(6)*(2+xB(4)/xB(5))/6 ...
%     + panel2Length(5)*panel2Thickness(5)*(2+xB(14)/xB(5))/6 ...
%     + panel3Length(6)*panel3Thickness(6)*(2+xB(6)/xB(5))/6;
% 
% A_effectiveD(6) = A_stringer + panel3Length(6)*panel3Thickness(6)*(2+xB(5)/xB(6))/6 + panel3Length(7)*panel3Thickness(7)*(2+xB(7)/xB(6))/6;
% A_effectiveD(7) = A_stringer + panel3Length(7)*panel3Thickness(7)*(2+xB(6)/xB(7))/6 + panel3Length(8)*panel3Thickness(8)*(2+xB(8)/xB(7))/6;
% A_effectiveD(8) = A_stringer + panel3Length(8)*panel3Thickness(8)*(2+xB(7)/xB(8))/6 + panel3Length(9)*panel3Thickness(9)*(2+xB(9)/xB(8))/6;
% 
% A_effectiveD(9) = A_spar_cap + panel3Length(9)*panel3Thickness(9)*(2+xB(8)/xB(9))/6 ...
%     + panel3Length(10)*panel3Thickness(10)*(2+xB(10)/xB(9))/6;

%% X CENTROID

B = [A_effective, fliplr(A_effective)]; % top half is equivlant to bottom half
% BD = [A_effectiveD, fliplr(A_effectiveD)];

xBCentroid = sum(xB.*B)/sum(B);
xBShifted = xB - xBCentroid; % shifts the x locations of the booms relative to the centroid 
%% BASIC SHEAR FLOW

Sy = SF; % vertical load applied
Sx = SF_D; % horizontal load applied
[Ixx, Iyy, Ixy] = MoACalc(B, yB, xBShifted);

qB1 = basicSFBoxCalc( B, yB, panel1Boom, Sy, Ixx, xBShifted, Sx, Iyy, Ixy);
qB2 = basicSFBoxCalc( B, yB, panel2Boom, Sy, Ixx, xBShifted, Sx, Iyy, Ixy);
qB3 = basicSFBoxCalc( B, yB, panel3Boom, Sy, Ixx, xBShifted, Sx, Iyy, Ixy);

%% RATE OF TWIST EQNS

% find common booms and panel indexes relating to the booms
[commonB1_2, panel1CommonIdx, panel2_1CommonIdx] = commonBooms(panel1Boom, panel2Boom);
[commonB3_2, panel3CommonIdx, panel2_3CommonIdx] = commonBooms(panel3Boom, panel2Boom);

% correction for multiple entries into shared panels
qB1(panel1CommonIdx) = qB1(panel1CommonIdx) + qB2(panel2_1CommonIdx-1);
qB2(panel2_1CommonIdx) = qB2(panel2_1CommonIdx) + qB1(panel1CommonIdx-1);
qB2(panel2_3CommonIdx) = qB2(panel2_3CommonIdx) + qB3(panel3CommonIdx-1);
qB3(panel3CommonIdx) = qB3(panel3CommonIdx) + qB2(panel2_3CommonIdx-1);

% rate of twist equation (Q11*qs01 + Q12*qs02 = Q13)
% cells 1 and 2
Q11 = 1/(2*A1*G)*sum(panel1Delta) + 1/(2*A2*G)*panel2Delta(panel2_1CommonIdx);
Q12 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) - 1/(2*A2*G)*sum(panel2Delta);
Q13 = 1/(2*A2*G)*panel2Delta(panel2_3CommonIdx);
S1 = - (1/(2*A1*G))*sum(qB1.*panel1Delta) + (1/(2*A2*G))*sum(qB2.*panel2Delta);

Q12 = Q12/Q11;
Q13 = Q13/Q11;
S1 = S1/Q11;
Q11 = Q11/Q11;

% cells 2 and 3
Q21 = -1/(2*A2*G)*panel2Delta(panel2_1CommonIdx);
Q22 = 1/(2*A2*G)*sum(panel2Delta) + 1/(2*A3*G)*panel2Delta(panel2_3CommonIdx);
Q23 = -1/(2*A2*G)*panel3Delta(panel3CommonIdx) - 1/(2*A3*G)*sum(panel3Delta);
S2 = - (1/(2*A2*G))*sum(qB2.*panel2Delta) + (1/(2*A3*G))*sum(qB3.*panel3Delta); 

% cells 1 and 3
% Q21 = 1/(2*A1*G)*sum(panel1Delta);
% Q22 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) + 1/(2*A3*G)*panel3Delta(panel3CommonIdx);
% Q23 = -1/(2*A3*G)*sum(panel3Delta);
% S2 = 1/(2*A3*G)*sum(panel3Delta.*qB3) - 1/(2*A1*G)*sum(panel1Delta.*qB1);

Q22 = Q22/Q21;
Q23 = Q23/Q21;
S2 = S2/Q21;
Q21 = Q21/Q21;
%% MOMENT EQUILIBRIUM
% moments are taken about the quarter chord centreline where both shear forces act through

% moment equilibrium (Q21*qs01 + Q22*qs02 = Q23)
Q31 = 2*A1;
Q32 = 2*A2;
Q33 = 2*A3;

% need to correct for basic shear flow signs inconsistency due to sign
% convention rotating with the cell (opposite on top to bottom - same moment)
S3 = -sum(qB1.*panel1Length.*panel1MomentArm) -sum(abs((qB2(1:end-1)).*panel2Length(1:end-1).*panel2MomentArm(1:end-1))) - sum(-abs(qB3(1:end-1).*panel3Length(1:end-1).*panel3MomentArm(1:end-1))) - qB3(end)*panel3Length(end)*panel3MomentArm(end);
% S3 = S3 - qB2(end)*panel2Length(end)*panel2MomentArm(end); % doubling up on stiffner 1-18

Q32 = Q32/Q31;
Q33 = Q33/Q31;
S3 = S3/Q31;
Q31 = Q31/Q31;


%% CONSTANT SHEAR FLOW
% solving the system of linear equations to provide the constant shear flow

Q = [Q11, Q12, Q13; Q21, Q22, Q23; Q31, Q32, Q33];
S = [S1; S2; S3];

qs0 = Q\S;

%% TOTAL SHEAR FLOW
% total shear flow is sum of basic and constant shear flow

qTotal1 = qB1 + qs0(1);
qTotal1(panel1CommonIdx) = qTotal1(panel1CommonIdx) - qs0(2);

qTotal2 = qB2 - qs0(2);
qTotal2(panel2_1CommonIdx) = qTotal2(panel2_1CommonIdx) + qs0(1);
qTotal2(panel2_3CommonIdx) = qTotal2(panel2_3CommonIdx) + qs0(3);

qTotal3 = qB3 + qs0(3);
qTotal3(panel3CommonIdx) = qTotal3(panel3CommonIdx) - qs0(2);

%% RATE OF TWIST EQUATIONS

dThetadZ1 = 1/(2*A1*G)*(sum(panel1Delta)*qs0(1) - panel1Delta(panel1CommonIdx)*qs0(2) + sum(qB1.*panel1Delta));
dThetadZ2 = 1/(2*A2*G)*(sum(panel2Delta)*qs0(2) - panel2Delta(panel2_1CommonIdx)*qs0(1) - panel2Delta(panel2_3CommonIdx)*qs0(3) + sum(qB2.*panel2Delta));
dThetadZ3 =  1/(2*A3*G)*(sum(panel3Delta)*qs0(3) - panel3Delta(panel3CommonIdx)*qs0(2) + sum(qB3.*panel3Delta));

%% VALIDATE RESULTS
% checks that result in common panels between cells are consistent
% ensures there are no major wrong equations
% does NOT ensure that shear flow values are correct (i.e. constant shear flow can change and not affect this as all cells are equally affected)

assert(abs((dThetadZ1 - dThetadZ2)/dThetadZ1) < PRECISE_RELATIVE_ERROR, 'Difference between rate of twist for cells 1 and 2 for shear');
assert(abs((dThetadZ2 - dThetadZ3)/dThetadZ2) < PRECISE_RELATIVE_ERROR, 'Difference between rate of twist for cells 2 and 3 for shear');
assert(abs((dThetadZ1 - dThetadZ3)/dThetadZ3) < PRECISE_RELATIVE_ERROR, 'Difference between rate of twist for cells 1 and 3 for shear');

debugTotal1 = qTotal1(panel1CommonIdx) - qTotal2(panel2_1CommonIdx);
assert(abs(debugTotal1/qTotal1(panel1CommonIdx)) < PRECISE_RELATIVE_ERROR, 'Difference in cells 1 and 2 common values');

debugTotal2 = qTotal3(panel3CommonIdx) - qTotal2(panel2_3CommonIdx);
assert(abs(debugTotal2/qTotal3(panel3CommonIdx)) < PRECISE_RELATIVE_ERROR, 'Difference in cells 2 and 3 common values');

%% SHEAR STRESS IN PANEL
% converts shear flow of each section to a shear stress

shearTotal1 = qTotal1./panel1Thickness;
shearTotal2 = qTotal2./panel2Thickness;
shearTotal3 = qTotal3./panel3Thickness;

%% TORSION SHEAR IN PANEL

T11 = 1/(2*A1*G)*sum(panel1Delta) + 1/(2*A2*G)*panel2Delta(panel2_1CommonIdx);
T12 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) - 1/(2*A2*G)*sum(panel2Delta);
T13 = 1/(2*A2*G)*panel2Delta(panel2_3CommonIdx);
TS1 = 0;

T12 = T12/T11;
T13 = T13/T11;
TS1 = TS1/T11;
T11 = T11/T11;

% cells 2 and 3
T21 = -1/(2*A2*G)*panel2Delta(panel2_1CommonIdx);
T22 = 1/(2*A2*G)*sum(panel2Delta) + 1/(2*A3*G)*panel2Delta(panel2_3CommonIdx);
T23 = -1/(2*A2*G)*panel3Delta(panel3CommonIdx) - 1/(2*A3*G)*sum(panel3Delta);
TS2 = 0; 

% cells 1 and 3
% Q21 = 1/(2*A1*G)*sum(panel1Delta);
% Q22 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) + 1/(2*A3*G)*panel3Delta(panel3CommonIdx);
% Q23 = -1/(2*A3*G)*sum(panel3Delta);
% TS2 = 0;

T22 = T22/T21;
T23 = T23/T21;
TS2 = TS2/T21;
T21 = T21/T21;

CM = -0.6; % low angle of attack < 20
T = q_MAX*AREA_WING*CM; % torque applied to wing in lb.in using the pitching moment coefficient of 0.5
T31 = 2*A1;
T32 = 2*A2;
T33 = 2*A3;
TS3 = T;

QT = [T11, T12, T13; T21, T22, T23; T13, T32, T33];

ST = [TS1; TS2; TS3];

qs0Torsion = QT\ST;

shearTorsion1 = qs0Torsion(1)./panel1Thickness;
shearTorsion2 = qs0Torsion(2)./panel2Thickness;
shearTorsion3 = qs0Torsion(3)./panel3Thickness;

dThetadZ1_torsion = 1/(2*A1*G)*(sum(panel1Delta)*qs0Torsion(1) - panel1Delta(panel1CommonIdx)*qs0Torsion(2));
dThetadZ2_torsion = 1/(2*A2*G)*(sum(panel2Delta)*qs0Torsion(2) - panel2Delta(panel2_1CommonIdx)*qs0Torsion(1) - panel2Delta(panel2_3CommonIdx)*qs0Torsion(3));
dThetadZ3_torsion =  1/(2*A3*G)*(sum(panel3Delta)*qs0Torsion(3) - panel3Delta(panel3CommonIdx)*qs0Torsion(2));

assert(abs((dThetadZ1_torsion - dThetadZ2_torsion)/dThetadZ1_torsion) < PRECISE_RELATIVE_ERROR, 'Difference between rate of twist for cells 1 and 2 for torsion');
assert(abs((dThetadZ2_torsion - dThetadZ3_torsion)/dThetadZ2_torsion) < PRECISE_RELATIVE_ERROR, 'Difference between rate of twist for cells 2 and 3 for torsion');
assert(abs((dThetadZ1_torsion - dThetadZ3_torsion)/dThetadZ3_torsion) < PRECISE_RELATIVE_ERROR, 'Difference between rate of twist for cells 1 and 3 for torsion');

%% COMBINED ROTATION

dThetadZ1_total = dThetadZ1 + dThetadZ1_torsion
dThetadZ2_total = dThetadZ2 + dThetadZ2_torsion;
dThetadZ3_total = dThetadZ3 + dThetadZ3_torsion;

thetaTotalWing = dThetadZ1_total*SPAN/2

%% COMBINED SHEAR

shearCombined1 = shearTotal1 + shearTorsion1;
shearCombined2 = shearTotal2 + shearTorsion2;
shearCombined3 = shearTotal3 + shearTorsion3;

%% AXIAL STRESS IN BOOMS
% calculates axial stress in the idealised booms

normalTotalB = -BM*yB/Ixx;
normalTotalD = -BM_D*xBShifted/Iyy;

%% COMBINED AXIAL STRESS

normalTotal = normalTotalB + normalTotalD;

[maxComp, BoomMaxComp] = min(normalTotal);
idxMaxComp = find(panel3Boom == BoomMaxComp);

%% SAVING RESULTS


if(SAVE_DATA == 1)

    format short g

    panel1Format = strcat(num2str(panel1Boom'), '-', num2str(circshift(panel1Boom',[-1 1])));
    panel2Format = strcat(num2str(panel2Boom'), '-', num2str(circshift(panel2Boom',[-1 1])));
    panel3Format = strcat(num2str(panel3Boom'), '-', num2str(circshift(panel3Boom',[-1 1])));

    namesBoom = {'BoomNumber'; 'BoomArea'; 'xLocation'; 'yLocation'; 'BAxStress'; 'DAxStress'; 'TAxStress'};
    tableBoom = table([1:18]', round(B', 4, 'significant'), round(xB', 4, 'significant'), round(yB', 4, 'significant'), round(normalTotalB'/KILO, 4, 'significant'), round(normalTotalD'/KILO, 4, 'significant'), round(normalTotal'/KILO, 4, 'significant'), 'VariableNames', namesBoom)

    format short
    namesPanel = {'PanelNumber'; 'Length'; 'Thickness'; 'Delta'; 'LDShear'; 'TShear'; 'CShear'};
    tablePanel1 = table(panel1Format, round(panel1Length', 4, 'significant'), round(panel1Thickness', 4, 'significant'), round(panel1Delta', 4, 'significant'), round(shearTotal1'/KILO, 4, 'significant'), round(shearTorsion1'/KILO, 4, 'significant'), round(shearCombined1'/KILO, 4, 'significant'), 'VariableNames', namesPanel)
    tablePanel2 = table(panel2Format, round(panel2Length', 4, 'significant'), round(panel2Thickness', 4, 'significant'), round(panel2Delta', 4, 'significant'), round(shearTotal2'/KILO, 4, 'significant'), round(shearTorsion2'/KILO, 4, 'significant'), round(shearCombined2'/KILO, 4, 'significant'), 'VariableNames', namesPanel)
    tablePanel3 = table(panel3Format, round(panel3Length', 4, 'significant'), round(panel3Thickness', 4, 'significant'), round(panel3Delta', 4, 'significant'), round(shearTotal3'/KILO, 4, 'significant'), round(shearTorsion3'/KILO, 4, 'significant'), round(shearCombined3'/KILO, 4, 'significant'), 'VariableNames', namesPanel)

    checkSave = input('Do you want to save the results (will overwrite any current data with same name in active directory)?', 's');
    checkSave = upper(checkSave);

    if (isequal(checkSave, 'Y') || isequal(checkSave, 'YES'))
        writetable(tableBoom, 'tableBoom.csv');
        writetable(tablePanel1, 'tablePanel1.csv');
        writetable(tablePanel2, 'tablePanel2.csv');
        writetable(tablePanel3, 'tablePanel3.csv');
    end

end

%% FAILURE CRITERION

tau = max(abs(shearCombined2));

F_VM = sqrt(max(abs(normalTotal))^2 + 3*(max(abs(shearTotal2)))^2);
F_VM_incTorsion = sqrt(max(abs(normalTotal))^2 + 3*(max(abs(shearCombined2)))^2);
FoS = F_YIELD/F_VM
FoS_incTorsion = F_YIELD/F_VM_incTorsion
FoSUlt_incTorsion = F_ULT/F_VM_incTorsion

%% FASTENER VALUES

N = max(normalTotal) * A_stringer

F_BR_YIELD = 64e3; 
F_BR_ULT = 118e3;
F_SHEAR = 41e3;

D_bolt = N/(F_BR_ULT*t_skin)

D_bolt = 1.13

Anet = N/F_ULT;
W = Anet/t_skin + D_bolt

L_joint = N/(2*t_skin*F_SHEAR)
e = L_joint + D_bolt/2;

Ashr = (N/(2*L_joint*t_skin))/F_SHEAR;
D2 = sqrt(4*Ashr/pi)

%% BUCKLING

k = 4; % conservatively assume all panels are simply supported
k_curv = 7;
k_stringer = 0.43;
v = 0.33; % Poisson's ratio of Al-2024
E_C = 75e9 * PA_TO_PSI; % compressive elastic modulus of Al-2024
boltSpacing = W; % NEED REFERENCE FOR THIS - contols inter-fastener buckling
boltSpacing = 2.8;
c = 1.5; % constants for counter shear rivets
NUM_RIBS = 16;
Le = (SPAN/2)/(NUM_RIBS+1) % effective length of sheets in wing acting as columns - controls stiffener column buckling

% skin/stiffener buckling (same equation)
FskinPanel1 = k*pi^2*E_C/(12*(1-v^2))*((panel1Thickness./panel1Length).^2);
FskinPanel2 = k*pi^2*E_C/(12*(1-v^2))*((panel2Thickness./panel2Length).^2);
FskinPanel3 = k*pi^2*E_C/(12*(1-v^2))*((panel3Thickness./panel3Length).^2);

FskinPanel1_curv = k_curv*pi^2*E_C/(12*(1-v^2))*((panel1Thickness./panel1Length).^2);
FskinPanel2_curv = k_curv*pi^2*E_C/(12*(1-v^2))*((panel2Thickness./panel2Length).^2);
FskinPanel3_curv = k_curv*pi^2*E_C/(12*(1-v^2))*((panel3Thickness./panel3Length).^2);

Re = sqrt(((LENGTH_1 + STIFFNER_1/2)/2)^2)
Rt_ratio = Re/panel1Thickness(1)
Z = panel1Length(1)^2/(Re*panel1Thickness(1)) * sqrt(1-v^2)
k_1 = 500
FskinPanel1(1) = k_1*pi^2*E_C/(12*(1-v^2))*((panel1Thickness(1)./panel1Length(1)).^2)

FoSSkin = FskinPanel3(idxMaxComp-1)/-maxComp
FoSSkin_curv = FskinPanel3_curv(idxMaxComp-1)/-maxComp

% stringer buckling
FstringerBuckling = k_stringer*pi^2*E_C/(12*(1-v^2))*(t_stringer/STRINGER_WIDTH)^2;
FoSStringer = FstringerBuckling/-maxComp

% inter fastener buckling
Fir1 = c*pi^2*E_C/(12*(1-v^2))*(panel1Thickness./boltSpacing).^2;
Fir2 = c*pi^2*E_C/(12*(1-v^2))*(panel2Thickness./boltSpacing).^2;
Fir3 = c*pi^2*E_C/(12*(1-v^2))*(panel3Thickness./boltSpacing).^2;

FoSir = Fir3(idxMaxComp - 1)/-maxComp

% stiffener column buckling
A_skin = 30*t_skin*t_skin;
y_skin = STRINGER_DEPTH + t_skin/2;
I_skin_xx = (1/12)*A_skin*t_skin^2;
I_skin_yy = (1/12)*(30*t_skin)^3*t_skin;

A_vertical = STRINGER_DEPTH*t_stringer;
x_vertical = 0;
y_vertical = STRINGER_DEPTH/2;
I_vertical_xx = (1/12)*t_stringer*STRINGER_DEPTH^3;
I_vertical_yy = (1/12)*STRINGER_DEPTH*t_stringer^3;

A_horizontal = (STRINGER_WIDTH-t_stringer)*t_stringer;
x_horizontal = (STRINGER_WIDTH-t_stringer)/2;
y_horizontal = STRINGER_DEPTH - t_stringer/2;
I_horizontal_xx = (1/12)*(STRINGER_WIDTH-t_stringer)*t_stringer^3;
I_horizontal_yy = (1/12)*t_stringer*(STRINGER_WIDTH-t_stringer)^3;

A_flange = (STRINGER_WIDTH-t_stringer)*t_stringer;
y_flange = t_stringer/2;
x_flange = -(STRINGER_WIDTH-t_stringer)/2;
I_flange_xx = (1/12)*t_stringer*STRINGER_WIDTH^3;
I_flange_yy = (1/12)*STRINGER_WIDTH*t_stringer^3;

% NO thin walled approximation
x_stiffener = (A_vertical*x_vertical + A_horizontal*x_horizontal + A_flange*x_flange)/(A_vertical + A_horizontal + A_flange); 
y_stiffener = (A_vertical*y_vertical + A_horizontal*y_horizontal + A_skin*y_skin + A_flange*y_flange)/(A_vertical + A_horizontal + A_skin + A_flange); 
I_stiffener_xx_NA = I_skin_xx + I_vertical_xx + I_horizontal_xx + I_flange_xx;
I_stiffener_xx = I_stiffener_xx_NA + A_skin*(y_skin - y_stiffener)^2 + A_vertical*(y_vertical - y_stiffener)^2 + A_horizontal*(y_horizontal - y_stiffener)^2 + A_flange*(y_flange - y_stiffener)^2;
I_stiffener_yy_NA = I_skin_yy + I_vertical_yy + I_horizontal_yy + I_flange_yy;
I_stiffener_yy = I_stiffener_yy_NA + A_vertical*(x_vertical - x_stiffener)^2 + A_horizontal*(x_horizontal - x_stiffener)^2 + A_flange*(x_flange - x_stiffener)^2;
I_stiffener = min(I_stiffener_xx, I_stiffener_yy); % buckles around least 2MoA
Pcr = pi^2*E_C*I_stiffener/(Le^2);
Fcr = Pcr/(A_vertical + A_horizontal + A_skin + A_flange);

FoScr = Fcr/-maxComp

% crippling
ce_2Free = 0.295;
ce_1Free = 0.317;
ce_0Free = 0.339;

btRatio = (STRINGER_WIDTH + STRINGER_DEPTH)/(2*t_stringer);
Fcc_stringer_base = ce_1Free*sqrt(F_YIELD*E_C)/(btRatio^0.75);
Fcc_stringer_top = ce_2Free*sqrt(F_YIELD*E_C)/(btRatio^0.75);
%Fcc_sparCap = ce_0Free*sqrt(F_YIELD*E_C)/(btRatio^0.75);

Fcc_stringer = min([Fcc_stringer_base, Fcc_stringer_top]);

FoScripPure = Fcc_stringer/-maxComp

if (Fcc_stringer > F_YIELD)
    Fcc_stringer = F_YIELD;
end

FoScripAdjusted = Fcc_stringer/-maxComp

end
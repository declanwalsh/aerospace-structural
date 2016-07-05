% Function that calculates the constant shear flow due 
% Author: Declan Walsh
% Last Modified: 5/7/2016

% Assumes that moment is taken in line with both shear forces and they don't contribute to moment at location

% INPUT ARGUMENTS
% panelBoom = booms (index) belonging to each panel in signed rotation order from cut
% qb = basic shear flow in each panel (following the panelBoom numbering)
% A = area of each cell
% G = shear modulus
% panelLength = length of each panel
% panelMomentArm = moment distance of panel from point moments are taken about

% RETURN OUTPUT
% qs0 = constant shear flow (3x1 array with cells 1, 2 and 3 constant shear flow)

function [ qs0 ] = constantSFCalcTorsion( panel1Boom, panel2Boom, panel3Boom, panel1Delta, panel2Delta, panel3Delta, qB1, qB2, qB3, A1, A2, A3, G, T)

    [Q11, Q12, Q13, S1, Q21, Q22, Q23, S2] = rateOfTwistCoefficients( panel1Boom, panel2Boom, panel3Boom, panel1Delta, panel2Delta, panel3Delta, qB1, qB2, qB3, A1, A2, A3, G);

    % MOMENT EQUILIBRIUM
    % moment equilibrium (Q21*qs01 + Q22*qs02 = Q23)
    Q31 = 2*A1;
    Q32 = 2*A2;
    Q33 = 2*A3;

    S3 = T; % only moment resulting is from torque due to no loads and no basic shear flow
    
    % normalising the result to the first term
    Q32 = Q32/Q31;
    Q33 = Q33/Q31;
    S3 = S3/Q31;
    Q31 = Q31/Q31;

    % CONSTANT SHEAR FLOW
    % solving the system of linear equations to provide the constant shear flow

    Q = [Q11, Q12, Q13; Q21, Q22, Q23; Q31, Q32, Q33];
    S = [S1; S2; S3];

    qs0 = Q\S;

end


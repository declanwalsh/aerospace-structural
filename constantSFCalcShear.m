% Function that calculates the constant shear flow due to shear loads (both x and y directions) both acting at the moment location
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
% panelDelta = length divided by thickness for each panel

% RETURN OUTPUT
% qs0 = constant shear flow (3x1 array with cells 1, 2 and 3 constant shear flow)

function [ qs0] = constantSFCalcShear( panel1Boom, panel2Boom, panel3Boom, qB1, qB2, qB3, A1, A2, A3, G, panel1Length, panel1MomentArm, panel1Delta, panel2Length, panel2MomentArm, panel2Delta, panel3Length, panel3MomentArm, panel3Delta)
    
    [Q11, Q12, Q13, S1, Q21, Q22, Q23, S2] = rateOfTwistCoefficients( panel1Boom, panel2Boom, panel3Boom, panel1Delta, panel2Delta, panel3Delta, qB1, qB2, qB3, A1, A2, A3, G);

    % MOMENT EQUILIBRIUM
    % moments are taken in line with both shear forces (i.e. they contribute no moment)

    % moment equilibrium (Q21*qs01 + Q22*qs02 = Q23)
    Q31 = 2*A1;
    Q32 = 2*A2;
    Q33 = 2*A3;

    % need to correct for basic shear flow signs inconsistency due to sign
    % convention rotating with the cell (opposite on top to bottom - same moment)
    S3 = -sum(qB1.*panel1Length.*panel1MomentArm) -sum(abs((qB2(1:end-1)).*panel2Length(1:end-1).*panel2MomentArm(1:end-1))) - sum(-abs(qB3(1:end-1).*panel3Length(1:end-1).*panel3MomentArm(1:end-1))) - qB3(end)*panel3Length(end)*panel3MomentArm(end);
    % S3 = S3 - qB2(end)*panel2Length(end)*panel2MomentArm(end); % doubling up on stiffner 1-18

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


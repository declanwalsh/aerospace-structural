% Calculates the rate of twist coefficients for a 3 cell wing box (i.e. gives you 2 equations in the form of their matricies)
% Author: Declan Walsh
% Last Modified: 5/7/2016

% INPUT ARGUMENTS
% panelBoom = booms (index) belonging to each panel in signed rotation order from cut
% qb = basic shear flow in each panel (following the panelBoom numbering)
% A = area of each cell
% G = shear modulus

% RETURN OUTPUT
% Q = coefficients
% S = numerical solutions

function [ Q11, Q12, Q13, S1, Q21, Q22, Q23, S2 ] = rateOfTwistCoefficients( panel1Boom, panel2Boom, panel3Boom, panel1Delta, panel2Delta, panel3Delta, qB1, qB2, qB3, A1, A2, A3, G )

    % find common booms and panel indexes relating to the booms
    [~, panel1CommonIdx, panel2_1CommonIdx] = commonBooms(panel1Boom, panel2Boom);
    [~, panel3CommonIdx, panel2_3CommonIdx] = commonBooms(panel3Boom, panel2Boom);

    % rate of twist equation (Q11*qs01 + Q12*qs02 = Q13)
    % cells 1 and 2
    Q11 = 1/(2*A1*G)*sum(panel1Delta) + 1/(2*A2*G)*panel2Delta(panel2_1CommonIdx);
    Q12 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) - 1/(2*A2*G)*sum(panel2Delta);
    Q13 = 1/(2*A2*G)*panel2Delta(panel2_3CommonIdx);
    S1 = - (1/(2*A1*G))*sum(qB1.*panel1Delta) + (1/(2*A2*G))*sum(qB2.*panel2Delta);

    % normalise with respect to first term
    Q12 = Q12/Q11;
    Q13 = Q13/Q11;
    S1 = S1/Q11;
    Q11 = Q11/Q11;

    % cells 2 and 3
    Q21 = -1/(2*A2*G)*panel2Delta(panel2_1CommonIdx);
    Q22 = 1/(2*A2*G)*sum(panel2Delta) + 1/(2*A3*G)*panel2Delta(panel2_3CommonIdx);
    Q23 = -1/(2*A2*G)*panel3Delta(panel3CommonIdx) - 1/(2*A3*G)*sum(panel3Delta);
    S2 = - (1/(2*A2*G))*sum(qB2.*panel2Delta) + (1/(2*A3*G))*sum(qB3.*panel3Delta); 

    % cells 1 and 3 - should give equivalent result
    % Q21 = 1/(2*A1*G)*sum(panel1Delta);
    % Q22 = -1/(2*A1*G)*panel1Delta(panel1CommonIdx) + 1/(2*A3*G)*panel3Delta(panel3CommonIdx);
    % Q23 = -1/(2*A3*G)*sum(panel3Delta);
    % S2 = 1/(2*A3*G)*sum(panel3Delta.*qB3) - 1/(2*A1*G)*sum(panel1Delta.*qB1);

    % normalise with respect to first term
    Q22 = Q22/Q21;
    Q23 = Q23/Q21;
    S2 = S2/Q21;
    Q21 = Q21/Q21;

end


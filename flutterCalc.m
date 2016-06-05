% Flutter calculator
% calculate bending and twisting frequencies at range of airspeeds and calculates when flutter occurs
% Author: Declan Walsh
% Last Modified: 5/06/2016

% based off Garth's flutter calculator excel spreadsheet

% Inputs:
% c = chord length
% freqFundTwist, freqFundBend = fundamental natural frequencies of twisting and bending of the structure with no airspeed
% a = percentage past mid chord where elastic axis is located

% Outputs:
% vFlutter = velocity at which flutter occurs (frequencies intersect)
% vDivergence = velocity at which divergence occurs (frequencies go to 0)

% Improvements can be made to generalise it

function [ vFlutter, vDivergence ] = flutterCalc( c, freqFundTwist, freqFundBend,a )

    RELATIVE_ERROR = 1e-3;

    b = c/2; % half chord length
    e = 0.25 + a; % percentage from mid chord where CG located
    
    r = sqrt(0.388); % I_p/(m*b^2)
    sigma = freqFundBend/freqFundTwist; % ratio of fundmental frequencies
    mu = 76; % m/(rho*pi*b^2)
    x = e - a; % percentage of half chord between EA and CG
    
    % bounds of calculations
    V_MIN = 0;
    V_MAX = 5.6; % from Garth's pretty plots
    % V_MAX = r*sqrt(mu/(1+2*a)) % from divergence speed
    NUM_POINTS = 100;
    
    airspeed = linspace(V_MIN, V_MAX, NUM_POINTS);
    
    V = airspeed/(b*freqFundTwist);
    s = zeros(length(V), 4);
    
    aEq = r^2 - x^2;
    bEq = r^2*(1+sigma^2) - airspeed.^2*((2/mu)*(0.5+a)+(2/mu)*x);
    cEq = sigma^2*r^2 - airspeed.^2*sigma^2*(2/mu)*(0.5+a);
    
    detEq = bEq.^2 - 4*aEq*cEq;
    
    % solves the determinant case
    s1 = sqrt(-bEq/(2*aEq) - sqrt(detEq)/(2*aEq));
    s2 = sqrt(-bEq/(2*aEq) + sqrt(detEq)/(2*aEq));
   
    LA = imag(s1) < 0;
    s1 = s1.*-LA + s1.*~LA; % converts if the imaginary number is negative
    
    figure;
    plot(airspeed, imag(s1), airspeed, imag(s2))
    title('Imaginary Solution - Frequencies')
    
    figure;
    plot(airspeed, real(s1), airspeed, real(s2))
    title('Real Solution - Magnitudes')
    
    vFlutter = airspeed(find((imag(s1) - imag(s2)) < RELATIVE_ERROR, 1));
    vDivergence = airspeed(find(imag(s1) < RELATIVE_ERROR, 1)); % may not be found in not enough data points (returns empty in that case)
    
%     syms s
%     for i = 1:length(V)
%         matrix = [s^2 + sigma^2, s^2*x + 2*V(i)^2/mu; s^2*x, s^2*r^2 + r^2 - (2*V(i)^2/mu)*(1/2 + a)];
%         dMatrix = det(matrix);
%         s(i,:) = solve(dMatrix, s)
%     end
        
end


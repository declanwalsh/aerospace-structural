% function that takes a general plane stress state in equilibrium and plots MC, finds principal stresses
% Author: Declan Walsh
% Last Modified: 20/03/2016

% Input Arguments
% Sxx = normal stress in x direction on element
% Sxy = normal stress in y direction on element
% Sxy = shear stress on element

function [ principalStresses ] = MC( Sxx, Syy, Sxy )

theta = linspace(0, 2*pi, 100); % angles for drawing the circle

S3 = 0;

avSxy = (Sxx + Syy)/2;
Rxy = sqrt((Sxx - avSxy)^2 + Sxy^2);
S1 = avSxy + Rxy;
S2 = avSxy - Rxy;
SxyMax = Rxy;
x1 = avSxy + Rxy*cos(theta);
y1 = Rxy*sin(theta);

Rxz = (avSxy - Rxy)/2;
avSxz = Rxz;
x2 = avSxz + Rxz*cos(theta);
y2 = Rxz*sin(theta);

Ryz = (avSxy + Rxy)/2;
avSyz = Ryz;
x3 = avSyz + Ryz*cos(theta);
y3 = Ryz*sin(theta);

% checks to validate

% plot the resultant Mohr's Circle
figure;
hold on
axis equal
axis on

% measured points
plot(Sxx, Sxy, 'k*');
plot(Syy, -Sxy, 'k*');

% principal stresses
plot(S1, 0, 'r*');
plot(S2, 0, 'r*');
plot(0, 0, 'r*');

% Mohr's circles
plot(x1, y1);
plot(x2, y2);
plot(x3, y3);

line([0 0], ylim, 'Color', 'k');  %x-axis
line(xlim, [0 0], 'Color', 'k');  %y-axis

principalStresses = [S1, S2, S3];

end
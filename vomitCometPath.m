v = 150; % m/s (assumed to be constant)
g = 9.81; % acceleration due to gravity

theta = linspace(-pi/4, pi/4, 50); % angles the aircraft rotates about
R = (v^2)./(g*cos(theta)); % solving for zero apparent weight of aircraft radius (i.e. weightlessness)

figure
subplot(1,3,1)
plot(theta, R)
title('Vomit Comet Path');
ylabel('Radius of turn (m)');
xlabel('Angle of aircraft (radians)');
subplot(1,3,2)
polar(theta, R);
subplot(1,3,3)
plot(R.*cos(theta), R.*sin(theta));
title('Vomit Comet Path');
xlabel('Horizontal Displacement (m)');
ylabel('Vertical Displacment (m)');
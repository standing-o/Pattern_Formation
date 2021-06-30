% This code follows the 'Solve Predator-Prey Equations' documentation of mathworks
% classical prey-predator model
% yp(1) = (1 - alpha*y(2))*y(1)
% yp(2) = (-1 + beta*y(1))*y(2)

clear; clc; close all;

% 0<t<15, x(0)=y(0)=20
% lotka function : alpha = 0.01, beta = 0.02
t0 = 0;
tfinal = 15;
y0 = [20; 20];   
[t,y] = ode23(@lotka,[t0 tfinal],y0);

% plot 1 : Populations Over Time
plot(t, y)
title('Predator/Prey Populations Over Time')
xlabel('t')
ylabel('Population')
legend('Prey','Predators','Location','North')

% plot 2 : Comparison between two groups
plot(y(:,1), y(:,2))
title('Phase Plane Plot')
xlabel('Prey Population')
ylabel('Predator Population')
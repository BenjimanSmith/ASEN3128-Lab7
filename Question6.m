% This script answers question 6 by simulating the longitudional dynamics
% of the 747-100 based off of the 4x4 A matrix derrived in Question3.m
% Author: Benjiman Smith
% Collaborators: E. Owen, I. Quezada
% Date: 3/06/2020
%
clear all
close all
W = 6.366e5;
W = 4.45*W; %N
g = 9.81; % m/s
m = W/g;

%% Table E1, Case II
Altitude = 20000; % ft
Altitude = Altitude*0.3048; % m
M = 0.5;
V = 518; % ft/s
V = V*0.3048; % m/s
W = 6.366e5; % lbf
W = 4.45*W; %N
Ix = 1.82e7; % slugs ft^2
Ix = 1.3558179619*Ix; % kg m^2
Iy = 3.31e7; % slugs ft^2
Iy = 1.3558179619*Iy; % kg m^2
Iz = 4.97e7; % slugs ft^2
Iz = 1.3558179619*Iz; % kg m^2
Izx = 9.70e5; % slugs ft^2
Izx = 1.3558179619*Izx; % kg m^2
zeta = deg2rad(-6.8); %rad
Cd = .040;
CaseII = [Altitude, M, V, W, Ix, Iy, Iz, Izx, zeta, Cd]; % define Case II SI unit equivalent of table E.1
%% Table E2, Longitudional
%1st column
Xu = -4.883e1; % lb*s/ft
Xu = Xu*(4.45/0.3048); % N*s/m
Xw = 1.546e3; % lb*s/ft
Xw = Xw*(4.45/0.3048); % N*s/m
Xq = 0; % lb*s/rad
Xwdot = 0; % lb*s^2/ft
Xdele = 3.994e4; % lb/rad
Xdele = Xdele*4.45; % N/rad
% 2nd column 
Zu = -1.342e3; % lb*s/ft
Zu = Zu*(4.45/0.3048); % N*s/m
Zw = -8.561e3; % lb*s/ft
Zw = Zw*(4.45/0.3048); % N*s/m
Zq = -1.263e5; % lb*s/rad
Zq = Zq*4.45; % N*s/rad
Zwdot = 3.104e2; % lb*s^2/ft
Zwdot = Zwdot*(4.45/0.3048); % N*s^2/m
Zdele = -3.341e5; % lb/rad
Zdele = Zdele*4.45; % N/rad
% 3rd column 
Mu = 8.176e3; % lb*s
Mu = Mu*4.45; % N/s
Mw = -5.627e4; % lb*s
Mw = Mw*4.45; % N/s
Mq = -1.394e7; % lb*s*ft/rad
Mq = Mq*4.45*0.3048; %N*s*m/rad
Mwdot = -4.138e3; % lb*s^2
Mwdot = Mwdot*4.45;
Mdele = -3.608e7; % lb*ft/rad
Mdele = Mdele*4.45*0.3048; % Nm/rad;
Converted = [Xu, Zu, Mu; Xw, Zw, Mw; Xq, Zq, Mq; Xwdot, Zwdot, Mwdot; Xdele, Zdele, Mdele]; % Table E.3 converted to SI units
%% Body frame to Stability frame
Xudot = Xu*((cos(zeta)^2)) - ((Xw + Zu)*sin(zeta)*cos(zeta)) + Zw*(sin(zeta)^2); % B.12,6
Xwprime = Xw*((cos(zeta)^2)) + ((Xu - Zw)*sin(zeta)*cos(zeta)) - Zu*(sin(zeta)^2); % B.12,6
Xqdot = (Xq*cos(zeta)) - (Zq*sin(zeta)); % B.12,6
Xudotdot = Zwdot*((sin(zeta))^2); % B.12,6
Xwdotdot = -Zwdot*sin(zeta)*cos(zeta); % B.12,6

Zudot = Zu*((cos(zeta)^2)) - ((Zw - Xu)*sin(zeta)*cos(zeta)) - Xw*(sin(zeta)^2); % B.12,6
Zwprime = Zw*((cos(zeta)^2)) + ((Zu + Xw)*sin(zeta)*cos(zeta)) + Xu*(sin(zeta)^2); % B.12,6
Zqdot = (Zq*cos(zeta)) + (Xq*sin(zeta)); % B.12,6
Zudotdot = -Zwdot*sin(zeta)*cos(zeta); % B.12,6
Zwdotdot = Zwdot*((cos(zeta))^2); % B.12,6

Mudot = (Mu*cos(zeta))-(Mw*sin(zeta)); % B.12,6
Mwprime = (Mw*cos(zeta)) + (Mu*sin(zeta)); % B.12,6
Mqdot = Mq; % B.12,6
Mudotdot = -Mwdot*sin(zeta); % B.12,6
Mwdotdot = Mwdot*cos(zeta); % B.12,6
Xcolumn = [Xudot, Xwprime, Xqdot, Xwdotdot]'; % define column of x values
Zcolumn = [Zudot, Zwprime, Zqdot, Zwdotdot]'; % define column of z values
Mcolumn = [Mudot, Mwprime, Mqdot, Mwdotdot]'; % define column of m values
DimDeriv = [Xcolumn, Zcolumn, Mcolumn]; % Matrix of stability derivatives in the stability frame


A = [Xudot/m, Xwprime/m ,0, -g*cos(0);...
    Zudot/(m-Zwdotdot), Zwprime/(m-Zwdotdot),  (Zqdot+m*V)/(m-Zwdotdot), -m*g*sin(0)/(m-Zwdotdot);...
    (Mudot +(Mwdotdot*Zudot/(m-Zwdotdot)))/Iy, (Mwprime +(Mwdotdot*Zwprime/(m-Zwdotdot)))/Iy, (Mqdot +(Mwdotdot*(Zqdot +m*V)/(m-Zwdotdot)))/Iy, -Mwdotdot*m*g*sin(0)/(Iy*(m-Zwdotdot));...
    0, 0, 1, 0];

%% Question 6

conditions = [0,0,0,0]'; % ?u^E is 10 m/s
tspan = linspace(0,300); % time vector from 0 300 seconds  
[t1, X] = ode45(@(t, F)Specs3LB5LC(t, F, A), tspan, conditions); % nonlinear ODE

%sublotting of all 8 variables
figure()

sgtitle('Linearized Lognitudional dynamics at trimmed state', 'interpreter', 'latex');
subplot(2,2,1);
plot(t1, X(:,1),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta u (m/s)')
title('change in deviation in u over time')

subplot(2,2,2);
plot(t1, X(:,2),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta w (m/s)')
title('change in w over time')

subplot(2,2,3);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

subplot(2,2,4);
plot(t1, X(:,4),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta \theta (rad)')
title('change in deviation in elevation over time')



% options = odeset('Events', @StopFnct, 'RelTol', 1e-8); % stop function that ends ODE when a tolerance of 1e-8 is met
conditions = [10,0,0,0]'; % ?u^E is 10 m/s
tspan = linspace(0,300); % time vector from 0 300 seconds  
[t1, X] = ode45(@(t, F)Specs3LB5LC(t, F, A), tspan, conditions); % nonlinear ODE

%sublotting of all 8 variables
figure()

sgtitle('Linearized Lognitudional dynamics with initial $\Delta u^{E}$ = 10 $\frac{m}{s}$', 'interpreter', 'latex');
subplot(2,2,1);
plot(t1, X(:,1),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta u (m/s)')
title('change in deviation in u over time')

subplot(2,2,2);
plot(t1, X(:,2),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta w (m/s)')
title('change in w over time')

subplot(2,2,3);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

subplot(2,2,4);
plot(t1, X(:,4),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta \theta (rad)')
title('change in deviation in elevation over time')

conditions = [0,10,0,0]'; % ?w^E is 10 m/s
tspan = linspace(0,300); % time vector from 0 300 seconds  
[t1, X] = ode45(@(t, F)Specs3LB5LC(t, F, A), tspan, conditions); % nonlinear ODE


figure()

sgtitle('Linearized Lognitudional dynamics with initial $\Delta w^{E}$ = 10 $\frac{m}{s}$', 'interpreter', 'latex');
subplot(2,2,1);
plot(t1, X(:,1),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta u (m/s)')
title('change in deviation in u over time')

subplot(2,2,2);
plot(t1, X(:,2),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta w (m/s)')
title('change in w over time')

subplot(2,2,3);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

subplot(2,2,4);
plot(t1, X(:,4),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta \theta (rad)')
title('change in deviation in elevation over time')
conditions = [0,0,0.1,0]'; % q^E is 0.1 rad/s
tspan = linspace(0,300); % time vector from 0 300 seconds  
[t1, X] = ode45(@(t, F)Specs3LB5LC(t, F, A), tspan, conditions); % nonlinear ODE

figure()

sgtitle('Linearized Lognitudional dynamics with initial $\Delta q$ = 0.1 $\frac{rad}{s}$', 'interpreter', 'latex');
subplot(2,2,1);
plot(t1, X(:,1),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta u (m/s)')
title('change in deviation in u over time')

subplot(2,2,2);
plot(t1, X(:,2),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta w (m/s)')
title('change in w over time')

subplot(2,2,3);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

subplot(2,2,4);
plot(t1, X(:,4),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta \theta (rad)')
title('change in deviation in elevation over time')

conditions = [0,0,0,0.1]'; % q is 0.1 rad/s
tspan = linspace(0,300); % time vector from 0 300 seconds  
[t1, X] = ode45(@(t, F)Specs3LB5LC(t, F, A), tspan, conditions); % nonlinear ODE

figure()

sgtitle('Linearized Lognitudional dynamics with initial $\Delta \theta$ = 0.1 $rad$', 'interpreter', 'latex');
subplot(2,2,1);
plot(t1, X(:,1),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta u (m/s)')
title('change in deviation in u over time')

subplot(2,2,2);
plot(t1, X(:,2),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta w (m/s)')
title('change in w over time')

subplot(2,2,3);
plot(t1, X(:,3),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('q (rad/s)')
title('change in pitch rate over time')

subplot(2,2,4);
plot(t1, X(:,4),'linewidth', 2);
grid on
xlabel('time (s)')
ylabel('\Delta \theta (rad)')
title('change in deviation in elevation over time')
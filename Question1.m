% This script answers question 1 of assignment 7 by converting tables E1
% and E3 from the textbook into SI units
%   Author: Benjiman Smith
%   Collaborators: E. Owen, I. Quezada
%   Date: 3/6/2020
%
close all
clear all
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
%% Table E3, Longitudional
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



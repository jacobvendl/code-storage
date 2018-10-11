%Jacob Vendl
%ASEN 3128 Assignment 1
%Due Jan 29th

clear all
close all
clc

%4.572 meters

%Intial conditions
Vx=20;                   %horizontal velocity
Vy=0;                   %vertical velocity
Vz=20;                  %upward velocity
V = sqrt(Vx^2 + Vy^2);   %magnitude of velocity
g = 9.81;                %m/s^2
c_d = 0.54;               %drag
A = pi*(0.75/2)^2;       %m^2 area of cross section
m = 0.6;                %kg mass
x(1) = 0;                %intial x postion
y(1) = 0;                %inital y postion
rho = 1.225;               %kg/m^3
t(1) = 0;                %intial time
dt = 0.001;               %s set the intervals at which time will be evalutated
i = 1;                   %counter
windspeed = 0;           %windspeed [m/s]
F=0;                     %lifting force

figure(1); hold on

%set initial conditons
x(1)=0;
y(1)=0;
z(1)=0;
x1(1)=0;
y1(1)=0;
z1(1)=0;
%call the external function that does flight
[x,y,z] = calcTrajec(x,y,z,Vx,Vy,Vz,dt,t,i,g,rho,c_d,A,m,V,windspeed,F);
[x1,y1,z1] = calcTrajec(x1,y1,z1,Vx,Vy,Vz,dt,t,i,g,rho,0,A,m,V,windspeed,F);

figure(1); hold on
plot3(x,y,z)
plot3(x1,y1,z1)
xlabel('x')
ylabel('y')
zlabel('z')

figure(1)
title('Displacement of a Basketball, with Drag')
xlabel('Downrange Displacement [m]')
ylabel('Horizontal Displacement [m]')
zlabel('Vertical Displacement [m]')
axis equal
grid on
grid minor
axis([0 100 -1 1 0 5])
view(40,10)





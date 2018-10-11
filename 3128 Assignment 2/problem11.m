%Jacob Vendl
%ASEN 3128 Assignment 2, Problem 10

%model quadcopter moving at 5m/s and maintaining azimuth of 0 deg.

clear all
close all
clc

param.eta=1e-3;
param.zeta=3e-3;
param.alpha=2e-6;
param.beta=1e-6;
param.k=.0024;
param.m=.068;
param.Ix=6.8e-5;
param.Iy=9.2e-5;
param.Iz=1.35e-4;
param.d=.06;
param.g=9.80665;

%calculate trim thrust
F=param.m*param.g;

x_e=0;
y_e=0;
z_e=-1;
psi=180;
theta=20; 
phi=20; 
u=5;
v=5;
w=0;
p=4;
q=3;
r=0;

param.forcecontrol = F; 

y0=[x_e y_e z_e psi theta phi u v w p q r];
[time,quadvec]=ode45(@(t,y)quadcopter(t,y,param), [0.0 10.0], y0);

x=quadvec(:,1);
y=quadvec(:,2);
z=quadvec(:,3);

figure(1); hold on
plot3(x,y,z)
plot3(x(1),y(1),z(1),'g*')
plot3(x(end),y(end),z(end),'r*')
grid on
xlabel('North')
ylabel('East')
zlabel('Down')
title('Flightpath Visualized')
legend('Flightpath','Starting Point','Ending Point')
view(-50,26)

figure(2);hold on
plot(time, x,'LineWidth',2)
plot(time, y,'LineWidth',2)
plot(time, -z,'LineWidth',2)
xlabel('Time (s)')
ylabel('Distance (m)')
title('Components of the Flightpath over Time')
legend('North','East','Height from Zero')
grid on
grid minor

figure(3); hold on
plot(time, quadvec(:,10),'LineWidth',2)
plot(time, quadvec(:,11),'LineWidth',2)
plot(time, quadvec(:,12),'LineWidth',2)
xlabel('Time (s)')
ylabel('Rate of Change (rad/s)')
title('Yaw, Pitch, and Roll over Time')
legend('Yaw','Pitch','Roll')
grid on
grid minor

figure(4); hold on
plot(time, quadvec(:,4),'LineWidth',2)
plot(time, quadvec(:,5),'LineWidth',2)
plot(time, quadvec(:,6),'LineWidth',2)
xlabel('Time (s)')
ylabel('Angle (rad)')
title('Euler Angles over Time')
legend('Psi','Theta','Phi')
grid on
grid minor

figure(5); hold on
plot(time, quadvec(:,7),'LineWidth',2)
plot(time, quadvec(:,8),'LineWidth',2)
plot(time, quadvec(:,9),'LineWidth',2)
xlabel('Time (s)')
ylabel('Speed (m/s)')
title('Directional Speed over Time')
legend('u','v','w')
grid on
grid minor

data = load('RSdata_Drone01_1325.mat');
X = data.rt_estim.signals.values(:,1);
Y = data.rt_estim.signals.values(:,2);
Z = data.rt_estim.signals.values(:,3);
Psi = data.rt_estim.signals.values(:,4);
Theta = data.rt_estim.signals.values(:,5);
Phi = data.rt_estim.signals.values(:,6);
t = data.rt_sensor.time;

figure(6); hold on
plot3(X,Y,Z)
grid on
xlabel('North')
ylabel('East')
zlabel('Down')
title('Flightpath Visualized')
legend('Flightpath')
view(-50,26)

figure(7); hold on
plot(t,X,'LineWidth',2)
plot(t,Y,'LineWidth',2)
plot(t,Z,'LineWidth',2)
legend('x','y','z')
xlabel('Time [s]')
ylabel('Displacement [m]')
title('Displacement of Copter over Time')
grid on 
grid minor

figure(8); hold on
plot(t,Psi,'LineWidth',2)
plot(t,Theta,'LineWidth',2)
plot(t,Phi,'LineWidth',2)
legend('Azimuth','Elevation','Bank')
xlabel('Time [s]')
ylabel('Angle [rad]')
title('Euler Angles over Time')
grid on 
grid minor



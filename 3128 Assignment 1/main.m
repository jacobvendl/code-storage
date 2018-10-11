%PURPOSE: The purpose of this script is to set input parameters, make an
%ode45 call, and plot results of the bottle rocket trajectory
%INPUTS: N/A
%OUTPUTS: N/A
%ASSUMPTIONS: There are corresponding functions to call, input parameters
%are known
%UID: e3cd71fd21a6
%DATE CREATED: 10 November 2016
%DATE MODIFIED: 20 November 2016

%housekeeping
clear all
close all
clc

%define all initla parameters
global rhoA R cDrag pBot pAmb g aB
rhoA = 0.961; %kg/m^3  air density at sea level
R = 287;        %J/kg*°K Universal gas constant
mB = 0.05;       %kg      mass of the botte
g=9.81;         %m/s^2   gravity
aB = pi*(.03/2)^2;    %cross-sectional area of the bottle
pAmb = 12.03;   %psi      Atmospheric pressure
pAmb = 6894.76*pAmb;  %Pa Convert to Pascals

%SET THE FOUR MAIN INPUT PARAMETERS
vA = .001;          %m^3
cDrag = 0.5;        %drag coefficient, 0.3 to 0.5
pBot = 65;      %psi      Pressure in the bottle
pBot = 6894.76*pBot + pAmb;    %Pa  Convert to Pascals
theta0=45*(pi/180); %45°,but in radians

%set initial values
v0=0.0; 
x0=20;
y0=0;
z0=20;

%mass velocity theta force drag
y0=[mB theta0 v0 x0 y0 z0];
%call ode45 once to get the values of everything as a function of time
[time,rocketvec]=ode45('Rocket', [0.0 10.0], y0);

%read data from what ode45 returned
mass=rocketvec(:,1);
theta=rocketvec(:,2);
vel=rocketvec(:,3);
x=rocketvec(:,4);
y=rocketvec(:,6);
z=rocketvec(:,6);

%plot the height of the rocket versus the distance travelled
figure(1)
plot(x,z,'LineWidth',2.0);
hold on

%find the maximum of distance and height
maxX = max(max(x));
maxZ = max(max(z));

%print out the maximum distance and height values
fprintf('The rocket reached a maximum distance of %.2f \n',maxX);
fprintf('The rocket reached a maximum height of %.2f \n',maxZ);








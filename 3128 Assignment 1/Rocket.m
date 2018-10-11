%PURPOSE: The purpopse of this function is to use equations given in the
%lab document to calculate the trajectory of a bottle rocket. 
%INPUTS: accepts initial conditons for 7 variables: mass of rocket, launch
%angle, velocity, x position, z position, volume of air, mass of air
%OUTPUTS: returns a matrix of values for all 7 variables as they change
%with time
%ASSUMPTIONS: Initial values are possible to obtain and formatted
%correctly.
%UID: e3cd71fd21a6
%DATE CREATED: 10 November 2016
%DATE MODIFIED: 20 November 2016

function [dydt] = Rocket(t,var)

%bring in all the global variables from the main function
global rhoA cDrag g aB ;

%standard ode45 thing here to set all the input data as nicer variables
m = var(1);
theta = var(2);
Vel = var(3);
x = var(4);
y = var(5);
z = var(6);

%set the change in non-complex variables
dydt(1)=0;
dydt(2) = -g*cos(theta)/Vel;    %change in theta
dydt(4) = Vel*cos(theta);       %change in x
dydt(5) = 0;                    %change in y
dydt(6) = Vel*sin(theta);       %change in z

%CHECK TO MAKE SURE WE HAVEN'T CRASH-LANDED YET
if z<=0
    for i=1:5
        dydt(i)=0;
    end
end

F=0;

%use drag force equation
D = rhoA/2 * (Vel.^2)*cDrag*aB;

%set change in velocity value
dydt(3) = F/m - D/m - g*sin(theta);

%transpose the 7 returned values. This is just an ode45 formality
dydt=dydt';
end

%Jake Vendl
%hw4 problem 1

clear all
close all
clc

%define the given vectors in the problem
r0 = [1770;850.8;-85.68];
v0 = [-0.6911;1.229;0.8143];
mu = 4902.799;

%calculate angular momentum
h = cross(r0,v0);
fprintf('h = %.2f km^2/s\n',norm(h))

%calculate eccentricity
e = cross(v0,h)/mu - r0/norm(r0);
fprintf('e = %.4f \n',norm(e))

%inclination, does not need a sign check
i = acosd(h(3)/norm(h));
fprintf('i = %.4f deg \n',i)

%find the vector direction of the line of nodes
n = cross([0;0;1],h);

Omega = acosd(n(1)/norm(n));
%but we need to perform a sign check on this guy. Basically, we're going to
%check if the ascending node is in the first two quadrants or the last
%one. The value already given by acosd defaults to the 1st,2nd quadrants.
%So we'll check the y-value of the line of nodes to orient this orbit. 

if n(2) < 0
    Omega = Omega + 180; %if the y-value was negative, add 180 to w to put the 
    %ascending node in the 3rd/4th quadrant
end
fprintf('RAAN = %.4f deg \n',Omega)

w = acosd( dot(r0,e) / norm(n) / norm(e));
%this value also needs a sign check. You can check either e or h's third
%component, the z-component, to see if periapsis is above or below the x-y
%plane and adjust the argument of periapsis accordingly if periapsis occurs
%below the plane. 
if e(3) < 0 
    w = w + 180;
end
fprintf('argument of periapsis = %.4f deg \n',w)

%calculate the initial true anomaly from the conic equation
theta0 = acosd(1/norm(e) * (norm(h)^2 / mu / norm(r0) - 1));
%perform a check on theta star initial by dotting r0 with v0 and checking
%the sign of that quantity. The dot product of r0 and v0 indicates which
%half of the orbit the s/c is in (top or bottom), with a positive value
%meaning top half and traveling away from periapsis

if dot(r0,v0) < 0 
    theta0 = -theta0; %if dot(r0,v0) was negative, switch signs on theta0
end
fprintf('theta*0 = %.4f deg \n',theta0)

%now move on to f and g functions part of problem, starting at 1c
a = norm(h)^2 / mu / (1-norm(e)^2);
fprintf('a = %.4f km \n',a)

%change to eccentric anomaly for time calculations
E0 = 2*atan(sqrt((1-norm(e))/(1+norm(e)))*tand(theta0/2)); %in radians

%find the intiial time of the orbit 
t0 = (E0 - norm(e)*sin(E0)) / sqrt(mu/a^3);
fprintf('t0 = %.2fs \n',t0)

%find the time after 30 minutes
t1 = t0+1800;

%set a Newton's method solver convergence
tol =1e-6;
E1 = pi/2;

%declare g, which is effectively just the difference between the LHS and
%RHS of Kepler's equation
g = E1 - norm(e)*sin(E1) - sqrt(mu/a^3)*t1;

%start a while loop that kinda just goes forever. The break will come from
%within the loop
while(abs(g)>=0)
    %calculate g, the difference between LHS and RHS
    g = E1 - norm(e)*sin(E1) - sqrt(mu/a^3)*t1;
    
    %check if  g makes it within the tolerance. if it does, you're done!
    if abs(g) <= tol
        break
    end
    %if g wasn't small enough, calculate g_prime in order to proceed with
    %Newton's method
    g_prime = 1-norm(e)*cos(E1);
    
    %declare a new E, which is the old E minus the error in the calculation
    %divided by the derivative
    E_new = E1 - g/g_prime;
    
    %set E equal to the new E so that the process can repeat
    E1 = E_new;
end

theta1 = 2*atand(sqrt((1+norm(e))/(1-norm(e))) * tan(E1/2));
fprintf('theta*1 = %.4f deg \n',theta1)

dtheta = theta1 - theta0;
fprintf('delta theta* = %.4f deg \n',dtheta)

%define p, the semi-latus rectum
p = norm(h)^2 / mu;

r1 = p / (1+norm(e)*cosd(theta1));
fprintf('r1 = %.4f km \n',r1)

%okay now I'm actually at the point where I can use an f and g function
%method

%calculate f/g functions and their derivatives based on lecture slides
%definitions
f = 1 - r1/p * (1 - cosd(dtheta));
g = norm(r0)*r1 / sqrt(mu*p) * sind(dtheta);
fdot = sqrt(mu/p) * tand(dtheta/2)*((1-cosd(dtheta))/p - 1/r1 - 1/norm(r0));
gdot = 1 - norm(r0)/p * (1-cosd(dtheta));

%now plug in to find rf
rf = f*r0 + g*v0;
vf = fdot*r0 + gdot*v0;

%print off the new position and velocity vectors 
fprintf('rf = %.3fx + %.3fy + %.3fz [km] \n',rf(1),rf(2),rf(3))
fprintf('vf = %.3fx + %.3fy + %.3fz [km/s] \n',vf(1),vf(2),vf(3))

hf = cross(rf,vf);
fprintf('hf = %.2f km^2/s (%.4f percent diff) \n',norm(hf),100*(norm(hf)-norm(h))/norm(h))

ef = cross(vf,hf)/mu - rf/norm(rf);
fprintf('ef = %.4f (%.4f percent diff)\n',norm(ef),100*(norm(ef)-norm(e))/norm(e))

af = norm(hf)^2 / mu / (1-norm(e)^2);
fprintf('af = %.4f km (%.4f percent diff) \n',af,(af-a)/af*100)

%calculate period of the orbit to be used in groundtrack analysis
T = 2*pi*sqrt(a^3/mu);
fprintf('T s/c = %.2f s \n',T)

omega_moon = 2*pi/(27.322*24*60*60) * 180/pi
delta_lon = omega_moon * T
    
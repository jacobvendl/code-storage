%Jake Vendl
%homework 4, problem 2

clear all
close all
clc

%define the given vectors in the problem
r1 = [-720000;670000;310000];
v1 = [2.16;-3.36;0.62];
mu = 3.794e7;

%calculate angular momentum
h = cross(r1,v1);
fprintf('h = %.2f km^2/s\n',norm(h))

%calculate eccentricity
e = cross(v1,h)/mu - r1/norm(r1);
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

w = acosd( dot(r1,e) / norm(n) / norm(e));
%this value also needs a sign check. You can check either e or h's third
%component, the z-component, to see if periapsis is above or below the x-y
%plane and adjust the argument of periapsis accordingly if periapsis occurs
%below the plane. 
if e(3) < 0 
    w = w + 180;
end
fprintf('argument of periapsis = %.4f deg \n',w)

%calculate the initial true anomaly from the conic equation
theta0 = acosd(1/norm(e) * (norm(h)^2 / mu / norm(r1) - 1));
%perform a check on theta star initial by dotting r0 with v0 and checking
%the sign of that quantity. The dot product of r0 and v0 indicates which
%half of the orbit the s/c is in (top or bottom), with a positive value
%meaning top half and traveling away from periapsis

if dot(r1,v1) < 0 
    theta0 = -theta0; %if dot(r0,v0) was negative, switch signs on theta0
end
fprintf('theta*1 = %.4f deg \n',theta0)

%calculate a, the semi-major axis
a = norm(h)^2 / mu / (1-norm(e)^2);
fprintf('a = %.4f km \n',a)

%change to eccentric anomaly for time calculations
E1 = 2*atan(sqrt((1-norm(e))/(1+norm(e)))*tand(theta0/2)); %in radians

%find the intiial time of the orbit 
t1 = (E1 - norm(e)*sin(E1)) / sqrt(mu/a^3);
fprintf('t1 = %.2fs \n',t1)

%find the time after 1 day
t2 = t1 + 3600*24;
fprintf('t2 = %.2fs \n',t2)

%set a Newton's method solver convergence
tol =1e-6;
E2 = pi/2;

%declare g, which is effectively just the difference between the LHS and
%RHS of Kepler's equation
g = E2 - norm(e)*sin(E2) - sqrt(mu/a^3)*t2;

%start a while loop that kinda just goes forever. The break will come from
%within the loop
while(abs(g)>=0)
    %calculate g, the difference between LHS and RHS
    g = E2 - norm(e)*sin(E2) - sqrt(mu/a^3)*t2;
    
    %check if  g makes it within the tolerance. if it does, you're done!
    if abs(g) <= tol
        break
    end
    %if g wasn't small enough, calculate g_prime in order to proceed with
    %Newton's method
    g_prime = 1-norm(e)*cos(E2);
    
    %declare a new E, which is the old E minus the error in the calculation
    %divided by the derivative
    E_new = E2 - g/g_prime;
    
    %set E equal to the new E so that the process can repeat
    E2 = E_new;
end

theta1 = 2*atand(sqrt((1+norm(e))/(1-norm(e))) * tan(E2/2));
fprintf('theta*2 = %.4f deg \n',theta1)

dtheta = theta1 - theta0;
fprintf('delta theta* = %.4f deg \n',dtheta)

%define p, the semi-latus rectum
p = norm(h)^2 / mu;

r2 = p / (1+norm(e)*cosd(theta1));
fprintf('r2 = %.4f km \n',r2)

%okay now I'm actually at the point where I can use an f and g function
%method

%calculate f/g functions and their derivatives based on lecture slides
%definitions
f = 1 - r2/p * (1 - cosd(dtheta));
g = norm(r1)*r2 / sqrt(mu*p) * sind(dtheta);
fdot = sqrt(mu/p) * tand(dtheta/2)*((1-cosd(dtheta))/p - 1/r2 - 1/norm(r1));
gdot = 1 - norm(r1)/p * (1-cosd(dtheta));

%now plug in to find rf
rf = f*r1 + g*v1;
vf = fdot*r1 + gdot*v1;

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


%PART C, Saturn impact
%find the true anomaly of the orbit when impact occurs
theta_sat = acosd(1/norm(e) * (p/60268-1));
theta_sat=-theta_sat; %manually choose the negative value

%change to eccentric anomaly for time calculations
E_sat = 2*atan(sqrt((1-norm(e))/(1+norm(e)))*tand(theta_sat/2)); %in radians

%find the time of impact by re-arranging Kepler's equation
t_impact = (E_sat-norm(e)*sin(E_sat)) / sqrt(mu/a^3);

%find the velocity of the impact using specific energy
v_impact = sqrt(2*mu/60268 - mu/a);



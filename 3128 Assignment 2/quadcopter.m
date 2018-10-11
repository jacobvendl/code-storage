function [dydt] = quadcopter(t,y,param)

f = param.forcecontrol/4*ones(1,4);

phi=y(4);
theta=y(5);
psi=y(6);
u=y(7);
v=y(8);
w=y(9);
p=y(10);
q=y(11);
r=y(12);

%define quadcopter-specific forces and moments
Xc = 0;
Yc = 0;
Zc = -(f(1)+f(2)+f(3)+f(4));
Lc = param.d/sqrt(2) * (f(2)+f(3)-f(1)-f(4));
Mc = param.d/sqrt(2) * (-f(1)-f(2)+f(3)+f(4));
Nc = param.k*(f(2)+f(4)-f(1)-f(3));

Xa = -param.eta*u^2*sign(u);
Ya = -param.eta*v^2*sign(v);
Za = -param.zeta*w^2*sign(w);
La = -param.alpha*p^2*sign(p);
Ma = -param.alpha*q^2*sign(q);
Na = -param.beta*r^2*sign(r);

%add together control and aerodynamic forces and moments
X=Xc+Xa;
Y=Yc+Ya;
Z=Zc+Za;
L=Lc+La;
M=Mc+Ma;
N=Nc+Na;

%use book equations 
xdot = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi)...
    -cos(phi)*sin(psi)) + w*(cos(phi)*sin(theta)*cos(psi)...
    +sin(phi)*sin(psi));
ydot = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi)...
    +cos(phi)*cos(psi)) + w*(cos(phi)*sin(theta)*sin(psi)...
    -sin(phi)*cos(psi));
zdot = -u*sin(theta) + v*sin(phi)*cos(theta) + w*cos(phi)*cos(theta);

phidot = p + tan(theta)*(q*sin(phi)+r*cos(phi));
thetadot = q*cos(phi) - r*sin(phi);
psidot = (q*sin(phi) + r*cos(phi))*sec(theta);

udot = X/param.m - param.g*sin(theta) - q*w + r*v;
vdot = Y/param.m + param.g*cos(theta)*sin(phi) - r*u + p*w;
wdot = Z/param.m + param.g*cos(theta)*cos(phi) - p*v + q*u;

pdot = L/param.Ix + ((param.Iy-param.Iz)/param.Ix)*q*r;
qdot = M/param.Iy + ((param.Iz-param.Ix)/param.Iy)*p*r;
rdot = N/param.Iz + ((param.Ix-param.Iy)/param.Iz)*p*q;

%delcare changes in time
dydt(1)=xdot;
dydt(2)=ydot;
dydt(3)=zdot;
dydt(4)=phidot;
dydt(5)=thetadot;
dydt(6)=psidot;
dydt(7)=udot;
dydt(8)=vdot;
dydt(9)=wdot;
dydt(10)=pdot;
dydt(11)=qdot;
dydt(12)=rdot;

%transpose because that's what ode45 makes you do
dydt=dydt';
end


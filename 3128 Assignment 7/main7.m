%Jake Vendl

clear all
close all
clc

%standard conversion factors
slug2kg = 14.5939;
ft2m=.3048;
lb2N = 4.4482;

%these are mainly pulled from Appendix E3,page 371
V = 871 * ft2m;
W = 6.366*10^5 * lb2N;
Ix = 1.82*10^7 * slug2kg * ft2m^2;
Iy = 3.31*10^7 * slug2kg * ft2m^2;
Iz = 4.97*10^7 * slug2kg * ft2m^2;
Izx = 9.70*10^5 * slug2kg * ft2m^2;
[~,~,~,rho]=atmoscoesa(40000*ft2m);
zeta = -2.4;
cD = 0.043;
S = 5500 * ft2m^2;
b = 195.68 * ft2m;
c_bar = 27.31 * ft2m;
h=0.25;

%defined according to Table 6.1 on page 166
Cxu = -.108;
Cxa = .2193;
Cxq = 0;
Cxa_dot = 0;
Czu = -.106;
Cza = -4.92;
Czq = -5.921;
Cza_dot = 5.896;
Cmu = .1043;
Cma = -1.023;
Cmq = -23.92;
Cma_dot = -6.314;
Cw0 = W / 0.5 / rho / V^2 / S; %defined on page 127
t_0 = 0; %given in Assignment 6 document
g=9.81;
m=W/g;

Xde = 1.544*10^4 * lb2N;
Zde = -3.677*10^5 * lb2N;
Mde = -4.038*10^7 * lb2N*ft2m;

%calculate dimensional derivatives, according to page 118
dim = [rho*V*S*Cw0*sind(t_0) + 0.5*rho*V*S*Cxu -rho*V*S*Cw0*cosd(t_0) + 0.5*rho*V*S*Czu 0.5*rho*V*c_bar*S*Cmu;
    0.5*rho*V*S*Cxa 0.5*rho*V*S*Cza 0.5*rho*V*c_bar*S*Cma;
    0.25*rho*V*c_bar*S*Cxq 0.25*rho*V*c_bar*S*Czq 0.25*rho*V*c_bar^2*S*Cmq;
    0.25*rho*c_bar*S*Cxa_dot 0.25*rho*c_bar*S*Cza_dot 0.25*rho*c_bar^2*S*Cma_dot];

%form A matrix, according to page 112
A = [dim(1,1) / m dim(2,1)/m 0 -g*cos(t_0);
    dim(1,2)/(m-dim(4,2)) dim(2,2)/(m-dim(4,2)) (dim(3,2)+m*V)/((m-dim(4,2))) -W*sin(t_0)/(m-dim(4,2));
    1/Iy*(dim(1,3)+dim(4,3)*dim(1,2)/(m-dim(4,2))) 1/Iy*(dim(2,3)+dim(4,3)*dim(2,2)/(m-dim(4,2))) 1/Iy*(dim(3,3)+dim(4,3)*(dim(3,2)+m*V)/(m-dim(4,2))) -dim(4,3)*W*sin(t_0)/Iy/(m-dim(4,2));
    0 0 1 0];

B = [Xde/m 0;
    Zde/(m-dim(4,2)) 0;
    Mde/Iy + dim(4,3)*Zde/(Iy*(m-dim(4,2))) 0;
    0 0];

[vec,val] = eig(A);

%the eigenvalues are the same for the first two and the last two, so just
%use (1,1) and (3,3) to check the norm
wn_1 = imag(val(1,1));
wn_2 = imag(val(3,3));

%define short mode as the one with the shortest period and highest
%oscillation frequency. In other words, if the norm of the first
%eigenvalues is smaller than the norm of the second two, then we know the
%first eigenvectors go with phugoid. This was defined in lecture.
if wn_1 < wn_2
    phu_1 = val(1,1);
    phu_2 = val(2,2);
    short_1 = val(3,3);
    short_2 = val(4,4);
    fprintf('Phugoid Eigenvalues: %.4f +/- %.4fi \n',real(val(1,1)),imag(val(1,1)))
    fprintf('Short Eigenvalues: %.4f +/- %.4fi \n',real(val(3,3)),imag(val(3,3)))
else
    phu_1 = val(4,4);
    phu_2 = val(3,3);
    short_1 = val(2,2);
    short_2 = val(1,1);
    fprintf('Phugoid Eigenvalues: %.4f +/- %.4fi \n',real(val(3,3)),imag(val(3,3)))
    fprintf('Short Eigenvalues: %.4f +/- %.4fi \n',real(val(1,1)),imag(val(1,1)))
end

wn_phu = norm(phu_1); %rad/sec
wn_short = norm(short_1); %rad/sec
fprintf('omega_n phugoid= %.4f \n',wn_phu)
fprintf('omega_n_short= %.4f \n',wn_short)

%defined using page 164 eqns
xi_phu = -real(phu_1)/wn_phu;
xi_short = -real(short_1)/wn_short;
fprintf('xi_phugoid= %.4f \n',xi_phu)
fprintf('xi_short= %.4f \n\n',xi_short)

%use approximation from class
A_approx = [dim(3,3)/Iy V*dim(2,3)/Iy; 1 0];
[vec_approx,val_approx]=eig(A_approx);
fprintf('****************************************\n')
fprintf('Approximation of A from Class \n')
disp(A_approx)

wn_short_approx=norm(val_approx(1,1));
xi_short_approx = -real(val_approx(1,1))/wn_short_approx;
fprintf('Using Approximation: \n')
fprintf('omega_n_short = %.4f \n',wn_short_approx)
fprintf('xi_short = %.4f \n\n',xi_short_approx)

%equation on page 163
T_phu = 2*pi / wn_phu;
fprintf('Phugoid Period (actual) = %.4fs \n',T_phu)

%equation on page 172
T_phu_Lancaster =  pi*sqrt(2)*V/g;
fprintf('Phugoid Period (Lancaster) = %.4fs \n\n',T_phu_Lancaster)

ks = linspace(1,3,201);
for i=1:length(ks)
    pitch_stiffness = -ks(i)*V*dim(2,3)/Iy;
    k2 = V*dim(2,3)*(1-ks(i))/Mde;
    wn = sqrt(pitch_stiffness);
    k1 = (2*xi_short*wn*Iy + dim(3,3))/Mde;
    
    short_approx = [(dim(3,3)-k1*Mde)/Iy (V*dim(2,3)-k2*Mde)/Iy;
        1 0];
    [vec1,val1]=eig(short_approx);
    evals1(i)=val1(1,1);
    evals2(i)=val1(2,2);
    
    k = [0 0 k1 k2;
        0 0 0 0];
    mat = A-B*k;
    [vec2,val2]=eig(mat);
    fullmodel1(i)=val2(1,1);
    fullmodel2(i)=val2(2,2);
    fullmodel3(i)=val2(3,3);
    fullmodel4(i)=val2(4,4);
end

figure; hold on
plot(evals1,'b','LineWidth',2)
plot(evals1(1),'g*')
plot(evals1(end),'r*')
plot(evals2,'b','LineWidth',2)
plot(evals2(1,1),'g*')
plot(evals2(end),'r*')
legend('First Eigenvalues','Start of 1st','End of 1st','Second Eigenvalues','Start of 2nd','End of 2nd')

title('Short Period Approximation Eigenvalues')
xlabel('Real')
ylabel('Imaginary')
grid on; grid minor

figure; hold on
plot(fullmodel1,'b','LineWidth',2)
plot(fullmodel3,'r','LineWidth',2)
plot(fullmodel2,'b','LineWidth',2)
plot(fullmodel4,'r','LineWidth',2)
legend('Short Period Evals','Phugoid Evals')
title('Closed Loop Eigenvalues')
xlabel('Real')
ylabel('Imaginary')
grid on; grid minor

%Now move on to the simulation. Recognize from lecture notes that the delta
%dot equations depend on the following variables, so we'll need to get them
%into the ode45 with the varialble param
clear ks k1 k2

%use values given in the assignment:
ks = 2;
pitch_stiffness = -ks*V*dim(2,3)/Iy;
k2 = V*dim(2,3)*(1-ks)/Mde;
wn = sqrt(pitch_stiffness);
if k2 ~= 0
    k1 = (2*xi_short*wn*Iy + dim(3,3))/Mde;
else
    k1=0;
end

param = [V t_0 m g Iy dim(1,1) dim(2,1) dim(3,1) dim(4,1)...
    dim(1,2) dim(2,2) dim(3,2) dim(4,2) dim(1,3)...
    dim(2,3) dim(3,3) dim(4,3) Xde Zde Mde k1 k2];

deltau=0;
deltaw=0;
deltaq=0;
deltatheta=0.1;
deltaxdot=0;
deltazdot=0;

initial_state=[deltau deltaw deltaq deltatheta deltaxdot deltazdot]';

%A(1,3:4)=0;
%A(2,3:4)=0;
%A(3,1:2)=0;
%A(4,1:2)=0;
%B(1,1)=0;
%B(2,1)=0;
%B(3,1)=0;
[t,planevec]=ode45(@(t,y) ode_plane(t,y,param,A,B),[0 500],initial_state);


%DISPLAY RESULTS
c = inv(vec)*[deltau deltaw deltaq deltatheta]';
fprintf('****************************************\n')
fprintf('With perturbation defined as: \n')
disp([deltau deltaw deltaq deltatheta]')
fprintf('The eigenvector coefficients are: \n')
disp(c)
if wn_1<wn_2
    if norm(c(1)) < norm(c(3))
        fprintf('This perturbtion mainly excites the short mode \n\n')
    else
        fprintf('This perturbtion mainly excites the phugoid mode \n\n')
    end
else
    if norm(c(1)) < norm(c(3))
        fprintf('This perturbtion mainly excites the phugoid mode \n\n')
    else
        fprintf('This perturbtion mainly excites the short mode \n\n')
    end
end

figure; hold on
suptitle(sprintf('Response of System Over Time with ks=%.0f',ks))
subplot(3,2,1)
plot(t,planevec(:,1))
xlabel('Time [s]')
ylabel('u velocity [m/s]')
grid on; grid minor

subplot(3,2,3)
plot(t,planevec(:,2))
xlabel('Time [s]')
ylabel('w velocity [m/s]')
grid on; grid minor

subplot(3,2,5)
plot(t,planevec(:,3))
xlabel('Time [s]')
ylabel('q [rad/s]')
grid on; grid minor

subplot(3,2,2)
plot(t,planevec(:,4))
xlabel('Time [s]')
ylabel('Thetadot [m/s^2]')
grid on; grid minor

subplot(3,2,4)
plot(t,planevec(:,5))
xlabel('Time [s]')
ylabel('X position [m]')
grid on; grid minor

subplot(3,2,6)
plot(t,planevec(:,6))
xlabel('Time [s]')
ylabel('Z position [m]')
grid on; grid minor


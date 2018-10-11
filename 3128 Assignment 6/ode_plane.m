function dydt = ode_plane(t,y,param)

V=param(1);
t_0=param(2);
m=param(3);
g=param(4);
Iy=param(5);
Xu=param(6); Xw=param(7); Xq=param(8); Xwdot=param(9);
Zu=param(10); Zw=param(11); Zq=param(12); Zwdot=param(13);
Mu=param(14); Mw=param(15); Mq=param(16); Mwdot=param(17);

deltau=y(1);
deltaw=y(2);
deltaq=y(3);
deltatheta=y(4);
deltax=y(5);
deltaz=y(6);

%page 109
udot=1/m*(Xu*deltau + Xw*deltaw) - g*cos(t_0)*deltatheta;
wdot=1/m*(Zu*deltau + Zw*deltaw + Zwdot*deltatheta + Zq*deltaq) - g*sin(t_0)*deltatheta + V*deltaq;
qdot=1/Iy*(Mu*deltau + Mw*deltaw + Mwdot*deltatheta + Mq*deltaq);
thetadot=deltaq;
xdot = deltau*cos(t_0) - V*deltatheta*sin(t_0) + deltaw*sin(t_0);
zdot = -deltau*sin(t_0) - V*deltatheta*cos(t_0) + deltaw*cos(t_0);

dydt=[udot wdot qdot thetadot xdot zdot]';
end


function dydt = ode_plane(t,y,param,A,B)

V=param(1);
t_0=param(2);
m=param(3);
g=param(4);
Iy=param(5);
Xu=param(6);    Xw=param(7);    Xq=param(8);    Xwdot=param(9);
Zu=param(10);   Zw=param(11);   Zq=param(12);   Zwdot=param(13);
Mu=param(14);   Mw=param(15);   Mq=param(16);   Mwdot=param(17);
Xde=param(18);  Zde=param(19);  Mde=param(20);
k1=param(21);   k2=param(22);

deltau=y(1);
deltaw=y(2);
deltaq=y(3);
deltatheta=y(4);
deltax=y(5);
deltaz=y(6);

deltaDeltadot = -k2*deltatheta - k1*deltaq;

%page 231, (7.7,6)
temp_vec = A*[deltau deltaw deltaq deltatheta]' + B(:,1)*deltaDeltadot;
udot = temp_vec(1);
wdot = temp_vec(2);
qdot = temp_vec(3);
thetadot = temp_vec(4);

xdot = deltau*cos(t_0) - V*deltatheta*sin(t_0) + deltaw*sin(t_0);
zdot = -deltau*sin(t_0) - V*deltatheta*cos(t_0) + deltaw*cos(t_0);

dydt=[udot wdot qdot thetadot xdot zdot]';
end


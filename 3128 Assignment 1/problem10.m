%Jacob Vendl
%ASEN 3128 Assignment 1
%Due Jan 29th

clear all
close all
clc

%Intial conditions
g = 9.81;                %m/s^2
c_d = 0.5;               %drag
rho_air = .961;               %kg/m^3
rho_he = 0.164;          %kg/m^3
dt = 0.01;               %s set the intervals at which time will be evalutated

for V_balloon = 100:100:10000 %run through different He displacements
    for wind = 0:1:20 %run through different windspeeds
        F=g*V_balloon*(rho_air-rho_he);     %calculate buoyant force
        A=pi*((3/4/pi)*V_balloon)^(2/3);    %calculate cross-sectional area for drag
        m = 0.50 + V_balloon*rho_he;                %kg mass
        
        %set initial conditions
        Vx=0;
        Vy=0;
        Vz=0;
        x(1) = 0;
        y(1) = 0;
        z(1) = 0;
        t(1) = 0;
        i=1;
        
        while min(t)< 60 %run for 60 seconds
            V=sqrt(Vx^2+Vy^2+Vz^2);
            t = t + dt;
            i = i + 1;
            
            %do this just like in the other functions
            if (Vx>wind)
                Dx=-0.5*rho_air*(abs(Vx-wind)^2)*c_d*A;
            else
                Dx=0.5*rho_air*(abs(Vx-wind)^2)*c_d*A;
            end
            Dy=-0.5*rho_air*(Vy^2)*c_d*A;

            Dz=-0.5*rho_air*(abs(Vz)^2)*c_d*A;
            
            dx = Dx/m;
            dy = Dy/m;
            dz = Dz/m -g + F/m;
            dVx = Vx + dx*dt;
            dVy = Vy + dy*dt;
            dVz = Vz + dz*dt;
            x(i) = x(i-1) + dVx*dt + 0.5*dx*dt^2;
            y(i) = y(i-1) + dVy*dt + 0.5*dy*dt^2;
            z(i) = z(i-1) + dVz*dt + 0.5*dz*dt^2;
            Vx = dVx;
            Vy = dVy;
            Vz = dVz;
        end
        %track the ascent angles
        angles(wind+1,V_balloon/100)=atand(z(end)/x(end));
        if sqrt(x(end)^2 + y(end)^2) < z(end) %then ascent angle is more than 45°
%             fprintf('PASSES! Windspeed = %.0fm/s \n',wind)
%             fprintf('He Volume = %.2fm^3 \n',V_balloon)
%             fprintf('Ascent Angle = %.1f° \n\n',atand(z(end)/x(end)))
            figure(1); hold on
            plot(wind,V_balloon,'go')
        else
%             fprintf('FAILS! Windspeed = %.0fm/s \n',wind)
%             fprintf('He Volume = %.2fm^3 \n',V_balloon)
%             fprintf('Ascent Angle = %.1f° \n\n',atand(z(end)/x(end)))
         end
        clear x y z t i
    end
end

%plot the ascent angles three different ways
figure(1); hold on
title('Ascent Angle Check, Comapring Windspeed and He Volume','fontSize',20)
xlabel('Windspeed [m/s]','fontSize',14)
ylabel('He Volume [m^3]','fontSize',14)

figure(2)
imagesc(rot90(angles))
title('Ascent Angle, Comparing He Volume and Windspeed','fontSize',20)
xlabel('Windspeed [m/s]')
xticks([1 3 5 7 9 11 13 15 17 19 21])
xticklabels({'0','2','4','6','8','10','12','14','16','18','20'})
ylabel('He Volume [m^3]')
yticks([1 11 21 31 41 51 61 71 81 91])
yticklabels({'10000','9000','8000','7000','6000','5000','4000','3000','2000','1000'})

figure(3)
contour(flipud(rot90(angles)),[30 35 40 45 50 55 60 65 70 75 80 85 90],'ShowText','On')
title('Ascent Angle, Comparing He Volume and Windspeed','fontSize',20)
xlabel('Windspeed [m/s]','fontSize',14)
ylabel('He Volume x100 [m^3]','fontSize',14)






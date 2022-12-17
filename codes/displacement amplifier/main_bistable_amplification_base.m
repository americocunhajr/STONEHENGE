%% Bistable EH with displacement amplification mechanism with base excitation
% main_bistable_amplification_base.m
% -------------------------------------------------------------------------
% Main file to perform numerical experiments of a bistable energy harvester
% with displacement amplifier forced by harmonic base excitation.
% 
% last update: 06/12/2022
%
% authors:
% João Pedro C. V. Norenberg (jp.norenberg@unesp.br)
% Americo Cunha Jr
% Piotr Wolszczak
% Grzegorz Litak
% -------------------------------------------------------------------------

clear
clc

% -------------------------------------------------------------------------
%% potential energy
xp = -.018:0.0005:.018;

k1 = -63.451;% -221.7168;
k3 = 634509;

PE = k1/2.*xp.^2 + k3/4.*xp.^4;

plot(xp,PE,'k','linewidth',2.5)
set(gca,'FontSize',16,'LineWidth',1.2,'xtick',[],'ytick',[],'Visible','off');
%% dimensionl parameters
m1  = 3*0.022;    % kg
m2  = 0.022;      % kg  
c1  = 0.1;        % Ns/m
c2  = 0.125;      % Ns/m
k1  = 500;        % N/m
k2  = -63.451;  % N/m
k2n = 634509;     % N/m^3
theta = -4.57e-3; % N/V
Rl    = 5e3;      % ohm  
Cp    = 4.3e-8;   % F
% -------------------------------------------------------------------------
%% base of excitation 
A = 3e-3;       % m
wext = 10;        % Hz  
% -------------------------------------------------------------------------
%% equation of motion
func = @(t,x)[x(2);
              (-c1*(x(2)-A*2*pi*wext*cos(2*pi*wext*t)) + c2*(x(4)-x(2)) ...
                    - k1*(x(1)-A*sin(2*pi*wext*t)) + k2*(x(3)-x(1)) +   ...
                    k2n*(x(3)-x(1))^3 + theta*x(5))/m1;
              x(4);
              (-c2*(x(4)-x(2)) -k2*(x(3)-x(1)) -k2n*(x(3)-x(1))^3 -     ...
                    theta*x(5))/m2;
              -1/(Rl*Cp)*x(5) - theta/Cp*(x(4)-x(2))];
% -------------------------------------------------------------------------
%% time integration
Fs = 4*128;                 % sampling frequency (Hz)
dt = inv(Fs);              % time period (s)
Nt = 100*Fs;                % number of samples
Td = (Nt-1)*dt;
tspan = linspace(0,Td,Nt); % time vector 
% -------------------------------------------------------------------------
%% integration: Runge-Kutta 
x0    = [0 0 0 0 0];
opts  = odeset('RelTol',1e-6,'AbsTol',1e-9);
[time,xout] = ode45(func,tspan,x0,opts);
% -------------------------------------------------------------------------
%% plotting

figure(1)
plot(time,xout(:,1))
title('displacement of amplification mechanism')

figure(2)
plot(time,xout(:,3))
title('displacement of bistable harvester')

figure(3)
plot(time,xout(:,5))
title('voltage')

figure(4)
plot(xout(3*end/4:end,3),xout(3*end/4:end,4))
title('displacement of bistable harvester')

%% Power mean
Power = xout(3*end/4:end,5).^2/Rl;
Mean_power = mean(Power);
RMS_power  = rms(Power);
rm    = m1/m2;
rk    = k1/k2;
fprintf('P_rms   P_mean   rm     rk \n')
fprintf('%.4f     %.4f   %.2f  %.2f \n',RMS_power,Mean_power,rm,abs(rk))
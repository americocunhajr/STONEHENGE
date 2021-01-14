% -----------------------------------------------------------------
%  main_pmeh_ivp.m
% ----------------------------------------------------------------- 
%  This is the main file for a program which response the time-
%  series for the piezo-magneto-elastic beam.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jpcvalese@gmail.com)
%        Americo Cunha (americo@ime.uerj.br)
%
%  last update: Oct 20, 2020
% -----------------------------------------------------------------

clc
clear
close all


% program header
% -----------------------------------------------------------
disp(' ---------------------------------------------------')
disp('            IVP Bistable Energy Harvester           ')
disp('                                                    ')
disp(' by                                                 ')
disp(' João Pedro Norenberg (Unesp)                       ')
disp(' jpcvalese@gmail.com                                ')
disp(' ---------------------------------------------------')
% -----------------------------------------------------------

Xpar.ksi    = 0.01;  % mechanical damping ratio
Xpar.chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
Xpar.lambda = 0.05;  % dimensionless time constant reciprocal
Xpar.kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
Xpar.f      = 0.083; % amplitude of excitation
Xpar.Omega  = 0.8;   % frequency of excitation
Xpar.beta   = 0.0;   % nonlinear electromechanical coupling term

% time interval integration 
tspan = 0:0.01:1000;

% initial condition
x0    = 1;
xdot0 = 0;
v0    = 0;
IC    = [x0 xdot0 v0];

% system response function 
[time,Y] = pmeh_eom(Xpar,IC,tspan);

% serial time response
Qdisp = Y(:,1);      % displacement-time of system
Qvelo = Y(:,2);      % velocity-time of system
Qvolt = Y(:,3);      % voltage-time of system

% post-processing
disp(' ---------------------------------------------------')
disp('            post-processing: plotingg               ')
disp(' ---------------------------------------------------')

% plot displament
figure(1)
plot(time,Qdisp,'b','LineWidth',1)
set(gca, 'FontName' , 'Arial' , 'FontSize' ,13 );
xlabel('     time     ', 'FontSize' , 15 , 'FontWeight' , 'bold');
ylabel(' displacement ', 'FontSize' , 15 , 'FontWeight' , 'bold');
grid

% plot voltage
figure(2)
plot(time,Qvolt,'b','LineWidth',1)
set(gca, 'FontName' , 'Arial' , 'FontSize' ,13 );
xlabel('    time  ', 'FontSize' , 15 , 'FontWeight' , 'bold');
ylabel('  voltage ', 'FontSize' , 15 , 'FontWeight' , 'bold');
ylim([-1.2 1.2])
grid

disp(' ---------------------------------------------------')
disp('               successfully finished                ')
disp(' ---------------------------------------------------')
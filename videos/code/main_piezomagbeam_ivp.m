% -----------------------------------------------------------------
%  main_piezomagbeam_ivp.m
%
%  This is the main file for a program that simulates the
%  nonlinear dynamics of a piezo-magneto-elastic beam,
%  which evolves according to the follwing system of
%  ordinary differential equations
%
%    d2x/dt2 + 2*ksi*dx/dt - 0.5*x*(1-x^2) - chi*v = f*cos(Omega*t)
%
%    dv/dt + lambda*v + kappa*dx/dt = 0
%
%        +
%
%    initial conditions,
%  
%  where
%  
%   x(t)   - dimensionless displacement of the beam tip
%   v(t)   - dimensionless voltage across the load resistance
%   t      - dimensionless time
%   ksi    - mechanical damping ratio
%   chi    - dimensionless piezoeletric coupling term (mechanical)
%   f      - dimensionless excitation amplitude
%   Omega  - dimensionless excitation frequency
%   lambda - dimensionless time constant reciprocal
%   kappa  - dimensionless piezoeletric coupling term (eletrical)
%  
%  Reference:
%  
%  A. Erturk, J. Hoffmann, and D. J. Inman
%  A piezomagnetoelastic structure for broadband vibration
%  energy harvesting
%  Applied Physics Letters
%  vol. 94 pp. 254102, 2009
% ----------------------------------------------------------------- 
%  programmers: 
%         Americo Cunha Jr (americo.cunhajr@gmail.com)
%         Joao Pedro Norenberg (jpcvalese@gmail.com)  
%
%  last update: Oct 19, 2020
% -----------------------------------------------------------------


clc
clear all
close all

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
f      = 0.147; % dimensionless excitation amplitude
Omega  = 0.8;   % dimensionless excitation frequency
beta   = 0.5;   % nonlinear term of electromechanical coupling 

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage

% integrate the initial value problem
% -----------------------------------------------------------
tic

% ODE solver optional parameters
opt = odeset('RelTol',1.0e-9,'AbsTol',1.0e-9);

% initial time of analysis
t0 = 0.0;

% final time of analysis
t1 = 300.0;

% time interval of analysis
tspan = t0:0.01:t1;

% initial conditions
IC = [x0; xdot0; v0];

% physical parameters vector
phys_param = [ksi chi f Omega lambda kappa beta];

% ODE solver Runge-Kutta45
[time,y] = ode45(@(t,y)piezomagbeam(t,y,phys_param),tspan,IC,opt);

% time series of dimensionless displacement
Qdisp = y(:,1);

% time series of dimensionless velocity
Qvelo = y(:,2);

% time series of dimensionless voltage
Qvolt = y(:,3);

% output power
[power,mean_power] = piezomagbeam_power(time,Qvolt,phys_param);

% number of dimensionless time steps
Ndt = length(time);

% number of steps for steady state
Nss = round(0.5*Ndt);
Nss1 = round(0.8*Ndt);
Nss2 = round(0.8*Ndt);


toc
% -----------------------------------------------------------
%%
% animate harvesting device dynamics
% ...........................................................

disp(' ')
disp(' --- animate harvesting device dynamics --- ');
disp(' --- inertial frame of reference --- ');
disp(' ');
disp('    ... ');
disp(' ');

Njump  =  1;
xmin   = -1;
xmax   =  1;
ymin   = -1;
ymax   =  1;


[fig00,F]  = plot_piezomagbeam_animation_inertial(time(1:Njump:end),...
                                     Qdisp(1:Njump:end),...
                                     Qvolt(1:Njump:end),...
                                    f,Omega,xmin,xmax,ymin,ymax,beta,IC);
disp(' ');
disp(' -- finished --');
disp(' ');
                        
%%
% animate harvesting device dynamics
% ...........................................................

disp(' ')
disp(' --- animate harvesting device dynamics --- ');
disp(' --- non-inertial frame of reference --- ');
disp(' ');
disp('    ... ');
disp(' ');

Njump  =  1;
xmin   = -1;
xmax   =  1;
ymin   = -1;
ymax   =  1;


[fig00,F]  = plot_piezomagbeam_animation_mobile(time(1:Njump:end),...
                                     Qdisp(1:Njump:end),...
                                     Qvolt(1:Njump:end),...
                                     f,Omega,xmin,xmax,ymin,ymax,beta,IC);

disp(' ');
disp(' -- finished --');
disp(' ');
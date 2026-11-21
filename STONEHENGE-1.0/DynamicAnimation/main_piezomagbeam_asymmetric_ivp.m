% -----------------------------------------------------------------
%  main_piezomagbeam_asymmetric_ivp.m
%
%  This is the main file for a program that simulates the
%  nonlinear dynamics of a piezo-magneto-elastic beam,
%  which evolves according to the follwing system of
%  ordinary differential equations
%
%    d2x/dt2 + 2*ksi*dx/dt - 0.5*x*(1+2*delta*x-x^2) - (1+beta*|x|)*chi*v 
%                                                         = f*cos(Omega*t)
%
%    dv/dt + lambda*v + (1+beta*|x|)*kappa*dx/dt = 0
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
%   beta   - nonlinear piezoelectric coupling
%   delta  - asymmetric coefficient of potential energy
%   phi    - bias angle
%  
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Jul 05, 2021
% -----------------------------------------------------------------
%% Processing
clc
clear all
close all

% defining parameter of the model
Xpar.ksi    = 0.01;  % mechanical damping ratio
Xpar.chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
Xpar.lambda = 0.05;  % dimensionless time constant reciprocal
Xpar.kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
Xpar.f      = 0.083; % dimensionless excitation amplitude
Xpar.Omega  = 0.8;   % dimensionless excitation frequency
Xpar.beta   = 0.0;   % nonlinear term of electromechanical coupling 
Xpar.delta  = 0.15;  % asymmetric coefficient of potential energy
Xpar.phi    = 10;    % bias angle

% initial time of analysis
t0 = 0.0;

% final time of analysis
t1 = 300.0;

% time interval of analysis
tspan = t0:0.01:t1;

% initial condition
x0    = 1.0;   % initial dimensionless displacement
xdot0 = 0.0;   % initial dimensionless velocity
v0    = 0.0;   % initial dimensionless voltage

IC    = [x0; xdot0; v0];

% ODE solver Runge-Kutta45
[time,y] = piezomagbeam_asymmetric(Xpar,IC,tspan);

% time series of dimensionless displacement
Qdisp = y(:,1);

% time series of dimensionless velocity
Qvelo = y(:,2);

% time series of dimensionless voltage
Qvolt = y(:,3);

% Sctruct input
series.time = time;
series.Disp = Qdisp;
series.Velo = Qvelo;
series.Volt = Qvolt;

% -----------------------------------------------------------
%% Plot: inertial frame of reference
% -----------------------------------------------------------

% animate harvesting device dynamics
% ...........................................................

disp(' ')
disp(' --- animate harvesting device dynamics --- ');
disp(' --- non-inertial frame of reference --- ');
disp(' ');
disp('    ... ');
disp(' ');

dim.xmin   = -1;
dim.xmax   =  1;
dim.ymin   = -1;
dim.ymax   =  1;

plot_piezomagbeam_animation_inertial_asymmetric(series, Xpar.f  , Xpar.Omega,...
                                       dim , Xpar.beta, IC, Xpar.phi,Xpar.delta);
disp(' ');
disp(' -- finished --');
disp(' ');
%% Plot: non-inertial frame of reference
% -----------------------------------------------------------

% animate harvesting device dynamics
% ...........................................................

disp(' ')
disp(' --- animate harvesting device dynamics --- ');
disp(' --- non-inertial frame of reference --- ');
disp(' ');
disp('    ... ');
disp(' ');

dim.xmin   = -1;
dim.xmax   =  1;
dim.ymin   = -1;
dim.ymax   =  1;
 
plot_piezomagbeam_animation_mobile_asymmetric(series,  f  , Omega,...
                                        dim , beta, IC, phi,delta);


disp(' ');
disp(' -- finished --');
disp(' ');

% -----------------------------------------------------------------
%  main_piezomagbeam_ivp.m
%
%  This is the main file for a program to do spectral analysis
%  of the nonlinear dynamics of a piezo-magneto-elastic beam,
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
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Sep 14, 2016
% -----------------------------------------------------------------


clc
clear all
close all


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' Piezo-Magneto-Elastic Beam Dynamics                ')
disp(' (spectral analysis)                                ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'piezomagbeam_fourier_f010_w08';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
f      = 0.10;  % dimensionless excitation amplitude
Omega  = 0.8;   % dimensionless excitation frequency

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage

% physical parameters vector
phys_param = [ksi chi f Omega lambda kappa x0 xdot0 v0];
% -----------------------------------------------------------


% integrate the initial value problem
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integration of the nonlinear dynamics --- ');
disp(' ');
disp('    ... ');
disp(' ');


% ODE solver optional parameters
opt = odeset('RelTol',1.0e-9,'AbsTol',1.0e-9);

% minimun frequency of the band
freq_min = 0.0;

% maximum frequency of the band
freq_max = 2.0;

% Nyquist sampling frequency
freq_samp = 2*freq_max;

% time step or sampling period
dt = 1/freq_samp;

% number of samples for signal copy
Nsamp = 1024;

% number of signal copies
Ncopies = 50;

% number of time steps
Ndt = Nsamp*(Ncopies+1);

% initial time of analysis
t0 = 0.0;

% final time of analysis
t1 = t0 + Ndt*dt;

% temporal mesh or instants of analysis
tspan = linspace(t0,t1,Ndt);
%tspan = [t0 t1];

% initial conditions
IC = [x0; xdot0; v0];

% ODE solver Runge-Kutta45
[time,y] = ode45(@(t,y)piezomagbeam(t,y,phys_param),tspan,IC,opt);

% time series of dimensionless displacement
Qdisp = y(:,1);

% time series of dimensionless velocity
Qvelo = y(:,2);

% time series of dimensionless voltage
Qvolt = y(:,3);

% number of dimensionless time steps
Ndt = length(time);

% number of steps for steady state
Nss1 = round(0.7*Ncopies)*Nsamp;
Nss2 = Ndt-Nsamp;
Nss3 = Nsamp+1;

toc
% -----------------------------------------------------------


% spectral analysis of the signals
% -----------------------------------------------------------

Nfft = 2^(nextpow2(Nsamp)+1);
%Nfft = Nsamp;

[psd_Qdisp,freq] = ...
 signal_psd(Qdisp(Nss3:Ndt),freq_max,freq_samp,Nfft,Nsamp,Ncopies);

[psd_Qvelo,freq] = ...
 signal_psd(Qvelo(Nss3:Ndt),freq_max,freq_samp,Nfft,Nsamp,Ncopies);

[psd_Qvolt,freq] = ...
 signal_psd(Qvolt(Nss3:Ndt),freq_max,freq_samp,Nfft,Nsamp,Ncopies);

psd_Qdisp_sg = sgolayfilt(psd_Qdisp,5,41);
psd_Qvelo_sg = sgolayfilt(psd_Qvelo,5,41);
psd_Qvolt_sg = sgolayfilt(psd_Qvolt,5,41);

% -----------------------------------------------------------


% save simulation results
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- saving simulation results --- ');
disp(' ');
disp('    ... ');
disp(' ');

save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------



% post processing
% -----------------------------------------------------------
post_piezomagbeam_fourier
% -----------------------------------------------------------

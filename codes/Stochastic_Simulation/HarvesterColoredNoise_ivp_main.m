% -----------------------------------------------------------------
%  HarvesterColoredNoise_ivp_main.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% -----------------------------------------------------------------
%  This script has a program to simulate the nonlinear dynamics
%  of a bistable piezo-magneto-elastic harvester which evolves 
%  according to the following system of differential equations
%
%  d2xdt2 + 2ksi dxdt - 0.5x(1-x^2) - chi v = fcos(Omega t) + noise
%
%  dvdt + lambda v + kappa dxdt = 0
%
%      +
%
%  initial conditions
%  
%  where
%  
%   x(t)   - dimensionless displacement of the beam tip
%   v(t)   - dimensionless voltage across the load resistance
%   t      - dimensionless time
%   ksi    - mechanical damping ratio
%   chi    - dimensionless piezoeletric coupling term (mechanical)
%   f      - dimensionless external excitation amplitude
%   Omega  - dimensionless external excitation frequency
%   noise  - dimensionless external excitation noise
%   lambda - dimensionless time constant reciprocal
%   kappa  - dimensionless piezoeletric coupling term (eletrical)
%  
%  References:
%  
%  [1] V. G. Lopes, J. V. L. L. Peterson, and A. Cunha Jr
%      The nonlinear dynamics of a bistable energy harvesting 
%      system with colored noise disturbances
%      Journal of Computational Interdisciplinary Sciences, 
%      vol. 10, pp. 125, 2019
%      https://hal.archives-ouvertes.fr/hal-02010224v2
%  
%  [2] A. Erturk, J. Hoffmann, and D. J. Inman
%      A piezomagnetoelastic structure for broadband vibration
%      energy harvesting
%      Applied Physics Letters
%      vol. 94 pp. 254102, 2009
%      https://doi.org/10.1063/1.3159815
% ----------------------------------------------------------------- 


clc
clear
close all


% program header
% -----------------------------------------------------------
disp('====================================================')
disp(' Bistable Piezo-Magneto-Elastic Energy Harvester    ')
disp(' by                                                 ')
disp(' Americo Cunha Jr                                   ')
disp('                                                    ')
disp(' Colored-noise driven stochastic nonlinear dynamics ')
disp('====================================================')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'HarvesterColoredNoise_corrlen';
%case_name = 'HarvesterColoredNoise_corrlen05';
%case_name = 'HarvesterColoredNoise_corrlen1';
%case_name = 'HarvesterColoredNoise_corrlen10';
%case_name = 'HarvesterColoredNoise_corrlen50';
%case_name = 'HarvesterColoredNoise_corrlen100';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% deterministic physical parameters
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
f      = 0.147; % dimensionless excitation amplitude
Omega  = 0.8;   % dimensionless excitation frequency

x0    = 1.0;  % dimensionless initial displacement
xdot0 = 0.0;  % dimensionless initial velocity
v0    = 0.0;  % dimensionless initial voltage

% physical parameters structure
param.ksi    = ksi;
param.chi    = chi;
param.lambda = lambda;
param.kappa  = kappa;
param.f      = f;
param.Omega  = Omega;
% -----------------------------------------------------------


% colored-noise external forcing
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- generating colored-noise external forcing --- ');
disp(' ');
disp('    ... ');
disp(' ');

% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% random process mean function
noise_mean = 0.0;

% standard deviation
noise_std = 0.05*f;

% correlation time (s)
tau_corr = 1.0;

% noise intensity
D = noise_std^2*tau_corr;

% initial time
t0 = 0.0;

% final time
t1 = t0 + 1000.0;

% number of samples for time step
Nsamp = 10;

% time step
dt = tau_corr/Nsamp;

% number of time steps
Ndt = round((t1-t0)/dt);

% domain upper limit
domain_upp = 0.5*(t1-t0);

% number of eigenpairs to be computed
Neig = 100;

% computing the autocorr function eigenpairs
[lambda_noise,phi_noise,time_noise] = ...
             fredholm_expcorr_eig(domain_upp,tau_corr,Neig,Ndt);

% temporal mesh vector
time_noise = time_noise + domain_upp;
         
% generate uncorrelated random variables
Y_noise = randn(1,Neig);

% multiply eigenvalues by noise intensity
lambda_noise = lambda_noise*D;

% representation energy level
energy_level = 0.95;

% number of eigenpairs used in KL expansion
[~,Nkl] = max(cumsum(lambda_noise)/sum(lambda_noise) >= energy_level);

% KL representation of random process (Ns x Ndt)
noise = noise_mean + ...
   phi_noise(:,1:Nkl)*sqrt(diag(lambda_noise(1:Nkl)))*Y_noise(:,1:Nkl)';

% physical parameters structure
param.noise      = noise;
param.time_noise = time_noise;

toc
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
%opt = odeset('RelTol',1.0e-9,'AbsTol',1.0e-9);

% time interval of analysis
%tspan = [t0 t1];
tspan = linspace(t0,t1,Ndt);

% initial conditions
IC = [x0; xdot0; v0];

% Runge-Kutta45 ODE solver
[time,Y] = ode45(@(t,y)HarvesterColoredNoise(t,y,param),tspan,IC);

% recover the solution at the desired instants
%Y = interp1(time,Y,time_noise);

% time series of dimensionless displacement
Qdisp = Y(:,1);

% time series of dimensionless velocity
Qvelo = Y(:,2);

% time series of dimensionless voltage
Qvolt = Y(:,3);

% output power
[power,power_mean] = HarvesterPower(time,Qvolt,lambda);

% external force
force_harmonic = f*cos(Omega*time_noise)';
force_total    = force_harmonic + noise;

% number of time steps
Ndt = length(time);

% number of steps for steady state
Nss = round(0.75*Ndt);


toc
% -----------------------------------------------------------

% statistical analysis
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- statistical analysis --- ');
disp(' ');
disp('    ... ');
disp(' ');

Nbins = 64;
Nksd  = 2*Nbins;

% histograms
[Qdisp_bins,Qdisp_freq] = randvar_pdf(Qdisp,Nbins);
[Qvelo_bins,Qvelo_freq] = randvar_pdf(Qvelo,Nbins);
[Qvolt_bins,Qvolt_freq] = randvar_pdf(Qvolt,Nbins);
[power_bins,power_freq] = randvar_pdf(power,Nbins);

% kernel density estimator
[Qdisp_ksd ,Qdisp_supp] = randvar_ksd(Qdisp,Nksd);
[Qvelo_ksd ,Qvelo_supp] = randvar_ksd(Qvelo,Nksd);
[Qvolt_ksd ,Qvolt_supp] = randvar_ksd(Qvolt,Nksd);
[power_ksd ,power_supp] = randvar_ksd(power,Nksd);

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
HarvesterColoredNoise_ivp_post
% -----------------------------------------------------------


% animate harvesting device dynamics
% ...........................................................
Njump = round(0.5*(1/dt));
vtitle = 'piezo-magneto-elastic energy harvesting device';
legend = ' rigid base frame of reference';
xmin   = -1;
xmax   =  1;
ymin   = -1;
ymax   =  1;
vname  = [num2str(case_name),'__animation'];
%fig00  = plot_piezomagbeam_animation( time(1:Njump:end),...
%                                     Qdisp(1:Njump:end),...
%                                vtitle,legend,xmin,xmax,ymin,ymax);
%close(fig00);
% ...........................................................

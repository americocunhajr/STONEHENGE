
% -----------------------------------------------------------------
%  main_piezomagbeam_resocurves.m
%
%  This is the main file for a program that estimates the
%  resonance curves of a piezo-magneto-elastic beam,
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
%  last update: Dec 8, 2016
% -----------------------------------------------------------------


clc
clear all
close all


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' Piezo-Magneto-Elastic Beam Dynamics                ')
disp(' (resonance curves)                                 ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
%case_name = 'piezomagbeam_resocurves_f010';
%case_name = 'piezomagbeam_resocurves_f011';
case_name = 'piezomagbeam_resocurves_f012';

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
f      = 0.12;  % dimensionless excitation amplitude
Omega  = 0.0;   % dimensionless excitation frequency

x0    = 1.0;  % initial dimensionless displacement
xdot0 = 0.0;  % initial dimensionless velocity
v0    = 0.0;  % initial dimensionless voltage


if strcmp(case_name,'piezomagbeam_resocurves_f010')
    
    % dimensionless excitation amplitude
    f = 0.10;
    
elseif strcmp(case_name,'piezomagbeam_resocurves_f011')
    
    % dimensionless excitation amplitude
    f = 0.11;
    
elseif strcmp(case_name,'piezomagbeam_resocurves_f012')
    
    % dimensionless excitation amplitude
    f = 0.12;
    
end
% -----------------------------------------------------------


% define the frequency window of analysis
% -----------------------------------------------------------
Omega_min = 0.1; % min. of Omega
Omega_max = 3.0; % max. of Omega

N_w = 5; % number of points

% vector of frequencies
Omega = linspace(Omega_min,Omega_max,N_w);
% -----------------------------------------------------------


% estimate the resonance curves
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- estimation of the resonance curves --- ');
disp(' ');
disp('    ... ');
disp(' ');

t0 = 0.0;    % initial dimensionless time
t1 = 5.0e3;  % final dimensionless time

% initial conditions
IC = [x0; xdot0; v0];

% ODE solver optional parameters
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);


% preallocate memory for amplitude vectors
Adisp = zeros(1,N_w);
Avelo = zeros(1,N_w);
Avolt = zeros(1,N_w);


for iw = 1:N_w
    
    % physical parameters vector
    phys_param = [ksi chi f Omega(iw) lambda kappa x0 xdot0 v0];
    
    % ODE solver Runge-Kutta45
    [time,y] = ode45(@(t,y)piezomagbeam(t,y,phys_param),[t0 t1],IC,opt);

    % number of dimensionless time steps
    Ndt = length(time);

    % number of steps for steady state
    Nss = round(0.8*Ndt);
    
    % maximum amplitude of dimensionless displacement
	Adisp(iw) = max(abs(y(Nss:Ndt,1)));
    
    % maximum amplitude of dimensionless velocity
	Avelo(iw) = max(abs(y(Nss:Ndt,2)));
    
    % maximum amplitude of dimensionless voltage
	Avolt(iw) = max(abs(y(Nss:Ndt,3)));

end

toc
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
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');



% resonance frequencies
% ...........................................................
[Adisp_pks,Adisp_Omega] = findpeaks(Adisp);
[Avolt_pks,Avolt_Omega] = findpeaks(Avolt);

disp(' displacement resonace frequencies')
disp(' ')
disp(Omega(Adisp_Omega));
disp(' ')

disp(' voltage resonace frequencies')
disp(' ')
disp(Omega(Avolt_Omega));
disp(' ')
% ...........................................................



% plot resonance curves
% ...........................................................
gtitle = ' resonance curves';
leg1   = ' displacement';
leg2   = ' voltage';
xlab   = ' excitation frequency';
ylab   = ' maximum amplitude';
xmin   = 0.0;
xmax   = Omega_max;
ymin   = 0.0;
ymax   = round(max(max(Adisp_pks),max(Avolt_pks)));
gname  = num2str(case_name);
flag   = 'eps';
fig1  = plot_resocurves(Omega,Adisp,...
                        Omega,Avolt,...
                        Omega(Adisp_Omega(2:end)),...
                        gtitle,leg1,leg2,...
                        xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ...........................................................


toc
% -----------------------------------------------------------

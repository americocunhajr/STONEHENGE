
% -----------------------------------------------------------------
%  main_piezomagbeam_opt_ds.m
%
%  This is the main file for a program which solve an optimization
%  problem that seaks to maximize the output power of a 
%  piezo-magneto-elastic beam, which evolves according to the 
%  follwing system of ordinary differential equations
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
clear
close all


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' Piezo-Magneto-Elastic Beam Dynamics                ')
disp(' (optimization problem)                             ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp('                                                    ')
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
f      = 0.0;   % dimensionless excitation amplitude
Omega  = 0.0;   % dimensionless excitation frequency

x0    = 1.0;  % dimensionless initial displacement
xdot0 = 0.0;  % dimensionless initial velocity
v0    = 0.0;  % dimensionless initial voltage
%---------------------------------------------------------------


% 0-1 test for chaos parameters
% -----------------------------------------------------------

% tolerance for 0-1 test
tol01 = 0.1;

% parameter c for 0-1 test
cmin = 0.0;
cmax = 2*pi;

% number of test repeats
Nc = 100;

% oversampling flag
OSflag = 0;
% -----------------------------------------------------------


% define ODE solver parameters
% -----------------------------------------------------------

% ODE solver optional parameters
%opt = odeset('RelTol',1.0e-9,'AbsTol',1.0e-6);

% inicial condition vector
IC = [x0; xdot0; v0];

% initial time of analysis
t0 = 0.0;

% final time of analysis
t1 = 2500.0;

% time interval of analysis
tspan = [t0 t1];

% number of dimensionless time steps
%Ndt = length(tspan);

% number of steps for steady state
%Nss = round(0.5*Ndt);

% number of jumps
%Njump = round((Ndt-Nss)/Ndt01);
% -----------------------------------------------------------



% solve optimization problem via direct search
%---------------------------------------------------------------
tic

disp(' '); 
disp(' --- optimization problem via direct search --- ');
disp(' ');
disp('    ... ');
disp(' ');


% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setDefaultStream(rng_stream); % Matlab 2009
%RandStream.setGlobalStream(rng_stream); % Matlab 2013

% number of points
Np1 = 256;
Np2 = 256;

% penalty parameter vector (K01 p1max p1min p2max p2min)
Hpenalty = 10.0;

% control parameter 1
%p1_min = 0.08; %camplitude f
%p1_max = 0.10; % amplitude f
p1_min = 0.04; % coupling chi
p1_max = 0.06; % coupling chi

% control parameter 2
%p2_min = 0.75; % frequency w
%p2_max = 0.85; % frequency w
p2_min = 0.4; % coupling kappa
p2_max = 0.6; % coupling kappa

% vector of p1 samples (excitation amplitude)
p1 = linspace(p1_min,p1_max,Np1);

% vector of p2 samples (excitation frequency)
p2 = linspace(p2_min,p2_max,Np2);

% preallocate memory for PerfFunc
% (remember: x = columns and y = lines)
S = zeros(Np2,Np1);

% preallocate memory for power function
% (remember: x = columns and y = lines)
power = zeros(Np2,Np1);

% preallocate memory for 0-1 classifier
% (remember: x = columns and y = lines)
K01 = zeros(Np2,Np1);

% initialize PerfFunc maximum/mininum value
S_max_overall = -Inf;
S_min_overall = +Inf;

% initialize power maximum/minimum value
power_max_overall = -Inf;
power_min_overall = +Inf;

% initialize K01 maximum/minimum value
K01_max_overall = -Inf;
K01_min_overall = +Inf;

% initialize power maximum/minimum points
p1_max_overall    = 0.0;
p2_max_overall    = 0.0;
p1_min_overall    = 0.0;
p2_min_overall    = 0.0;


% simulation information
%case_name  = ['piezomagbeam_opt_ds_',...
%              'penal_',num2str(Hpenalty),'_',...
%                'Np1_',num2str(Np1),'_',...
%                'Np2_',num2str(Np2)];

case_name  = ['piezomagbeam_optdesign_ds_',...
              'penal_',num2str(Hpenalty),'_',...
                'Np1_',num2str(Np1),'_',...
                'Np2_',num2str(Np2)];
            

% define data file name
file_name = [case_name,'.dat'];

% open data file
fileID = fopen(file_name,'w');

% define data format
formatSpec = '%.4f %.4f %+0.4E %0.4E %1.4f \n';

% display extreme values on screen
%fprintf('\n\nLocal values: f Omega PerfFunc power K01\n\n');

for np1 = 1:Np1
    
    % print loop indicador
    disp(['Np1 = ',num2str(np1)]);
    
    for np2 = 1:Np2
        
        % update parameters
        %f     = p1(np1);
        %Omega = p2(np2);
        chi   = p1(np1);
        kappa = p2(np2);
        
        % define physical paramters vector
        phys_param = [ksi chi f Omega lambda kappa];
                          
        % penalized performance function S(x)
        [S(np2,np1),power(np2,np1),K01(np2,np1)] = ...
          piezomagbeam_opt_PerfFunc(phys_param,tspan,IC,...
                                    cmin,cmax,Nc,tol01,OSflag,Hpenalty,...
                                    p1_min,p1_max,p2_min,p2_max);
        
        % update global maximum
        if S(np2,np1) > S_max_overall
            
              S_max_overall = S(np2,np1);
          power_max_overall = power(np2,np1);
            K01_max_overall = K01(np2,np1);
             p1_max_overall = p1(np1);
             p2_max_overall = p2(np2);
        end
        
        % update global minimum
        if S(np2,np1) < S_min_overall
            
             S_min_overall = S(np2,np1);
         power_min_overall = power(np2,np1);
           K01_min_overall = K01(np2,np1);
            p1_min_overall = p1(np1);
            p2_min_overall = p2(np2);
        end
        
        % print local values in data file
        fprintf(fileID,formatSpec,...
                 p1(np1),p2(np2),S(np2,np1),power(np2,np1),K01(np2,np1));
        
    end
end

% print global maximum
fprintf(fileID,formatSpec,...
        p1_max_overall,p2_max_overall,...
        S_max_overall,power_max_overall,K01_max_overall);
    
% print global minimum
fprintf(fileID,formatSpec,...
        p1_min_overall,p2_min_overall,...
        S_min_overall,power_min_overall,K01_min_overall);

% close data file
fclose(fileID);


% display global extreme values on screen
fprintf('\n\nGlobal maximum (f Omega PerfFunc power K01):\n\n');

fprintf(formatSpec,...
        p1_max_overall,p2_max_overall,...
        S_max_overall,power_max_overall,K01_max_overall);
    
fprintf('\n\nGlobal minimum (f Omega PerfFunc power K01):\n\n');

fprintf(formatSpec,...
        p1_min_overall,p2_min_overall,...
        S_min_overall,power_min_overall,K01_min_overall);

fprintf('\n');

disp('   ');
time_elapsed = toc
%---------------------------------------------------------------



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


% plot mean power countourmap
% ...........................................................
gtitle = ' performance function contour map';
%xlab   = ' excitation amplitude';
%ylab   = ' excitation frequency';
xlab   = ' mechanical coupling';
ylab   = ' electrical coupling';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
gname  = [case_name,'__perf_func'];
flag   = 'eps';
fig1   = graph_contour_pnt(p1,p2,S',...
                            p1_max_overall,p2_max_overall,...
                            gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ...........................................................


% plot mean power countourmap
% ...........................................................
gtitle = ' mean power contour map';
%xlab   = ' excitation amplitude';
%ylab   = ' excitation frequency';
xlab   = ' mechanical coupling';
ylab   = ' electrical coupling';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
gname  = [case_name,'__mean_power'];
flag   = 'eps';
%fig2   = graph_contour_pnt(p1,p2,power',...
%                            p1_max_overall,p2_max_overall,...
%                            gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig2   = graph_contour_pnt(p1,p2,power',...
                            p1_max_overall,p2_max_overall,...
                            gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

%close(fig2);
% ...........................................................


% plot 0-1 identification parameter
% ...........................................................
%xlab   = 'excitation amplitude';
%ylab   = 'excitation frequency';
xlab   = ' mechanical coupling';
ylab   = ' electrical coupling';
label0 = 'regular';
label1 = 'chaos';
gtitle = '0-1 test for chaos';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
zmin   = 0;
zmax   = 1;
gname  = [case_name,'__01test'];
flag   = 'eps';
fig3   = graph_binarymap(p1,p2,K01,...                          
                          p1_max_overall,p2_max_overall,...
                          gtitle,xlab,ylab,label0,label1,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
% fig3   = graph_binarymap(p1,p2,K01,...                          
%                             [],[],...
%                           gtitle,xlab,ylab,label0,label1,...
%                           xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3);
% ...........................................................


% plot 0-1 identification parameter ZOOM
% ...........................................................
%xlab   = 'excitation amplitude';
%ylab   = 'excitation frequency';
xlab   = ' mechanical coupling';
ylab   = ' electrical coupling';
label0 = 'regular';
label1 = 'chaos';
gtitle = '0-1 test for chaos';
xmin   = 0.099;
xmax   = 0.1;
ymin   = 0.77;
ymax   = 0.78;
gname  = [case_name,'__01test_zoom'];
flag   = 'eps';
fig4   = graph_binarymap(p1,p2,K01,...                          
                          p1_max_overall,p2_max_overall,...
                          gtitle,xlab,ylab,label0,label1,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig4);
% ...........................................................

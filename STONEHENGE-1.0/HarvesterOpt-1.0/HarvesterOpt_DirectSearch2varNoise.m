% -----------------------------------------------------------------
%  HarvesterOpt_DirectSearch2varNoise.m
% ----------------------------------------------------------------- 
%  This program solves an optimization problem that seaks to 
%  maximize the mean output power of a piezo-magneto-elastic
%  energy harvester system, defined by 
%  
%                  -- t0 + T
%                 |
%  power := 1/T   | lambda v(t)^2 dt
%                 |
%               -- t = t0
%  
%  in which the voltage come from the following dynamical system 
%
%  d2x/dt2 + 2*ksi*dx/dt - 0.5*x*(1-x^2) - chi*v = f*cos(Omega*t)
%  dv/dt   + lambda*v    + kappa*dx/dt           = 0
%          +
%  initial conditions,
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
%  The solution strategy is based on a direct search in a fine mesh,
%  where the optimum is found by exhaustion of the design points. 
%  
%  The design variables are: (f,Omega)
%  
%  Reference:
%  A. Cunha Jr
%  Enhancing the performance of a bistable energy harvesting 
%  device via the cross-entropy method (2020)
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunha@uerj.br
%
%  last update: October 23, 2020
% -----------------------------------------------------------------

clc
clear
close all


% program header
% -----------------------------------------------------------
disp(' ---------------------------------------------------')
disp(' Noisily Piezo-Magneto-Elastic Harvester Optimizer  ')
disp(' (direct search)                                    ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Cunha (UERJ)                               ')
disp(' americo.cunha@uerj.br                              ')
disp(' ---------------------------------------------------')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'HarvesterOpt_DS_2var_Noise';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define physical parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining model parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
%f      = 0.115; % dimensionless excitation amplitude
%Omega  = 0.8;   % dimensionless excitation frequency

x0    = 1.0;  % dimensionless initial displacement
xdot0 = 0.0;  % dimensionless initial velocity
v0    = 0.0;  % dimensionless initial voltage

toc
%---------------------------------------------------------------


% 0-1 test for chaos parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining 0-1 test for chaos parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% RNG seed (seed is fixed for reprodutibility)
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);

% tolerance for 0-1 test
tol01 = 0.1;

% parameter c for 0-1 test
cmin = 0.0;
cmax = 2*pi;

% number of test repeats
Nc = 100;

% oversampling flag
OSflag = 0;

toc
% -----------------------------------------------------------


% define ODE solver parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining ODE solver parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

% initial time of analysis
t0 = 0.0;

% final time of analysis
t1 = 2500.0;

% time interval of analysis
tspan = [t0 t1];

% inicial condition vector
IC = [x0; xdot0; v0];

% ODE solver optional parameters
%opt = odeset('RelTol',1.0e-9,'AbsTol',1.0e-6);

toc
% -----------------------------------------------------------



% solve optimization problem via direct search
%---------------------------------------------------------------
tic
disp(' '); 
disp(' --- optimization via direct search --- ');
disp(' ');
disp('    ... ');
disp(' ');

% number of design variables
Ndv = 2;

% number of points for design variables discretization
Np1 = 256;
Np2 = 256;

% number of points in the numerical grid
Ngrid = Np1*Np2;

% square root of Ngrid
sqrt_Ngrid = round(sqrt(Ngrid));

% penalty parameter
H = 10.0;

% limits for the design variables: (f,Omega)
X_min = [0.08 0.75];
X_max = [0.10 0.85];

% numerical grid for domain discretization
X1 = linspace(X_min(1),X_max(1),Np1);
X2 = linspace(X_min(2),X_max(2),Np2);

% preallocate memory for the ObjFunc
S = zeros(Ngrid,1);

% preallocate memory for power function
% (x = columns and y = lines)
power = zeros(Ngrid,1);

% preallocate memory for 0-1 classifier
% (x = columns and y = lines)
K01 = zeros(Ngrid,1);

% preallocate memory for the optimum point
X_opt = zeros(Ndv,1);

% initialize ObjFunc maximum
S_opt = -Inf;

% initialize power maximum
power_opt = -Inf;

% initialize K01 maximum
K01_opt = -Inf;

% hyper parameters
hyperparam.tspan  = tspan;
hyperparam.IC     = IC;
hyperparam.cmin   = cmin;
hyperparam.cmax   = cmax;
hyperparam.Nc     = Nc;
hyperparam.tol01  = tol01;
hyperparam.OSflag = OSflag;
hyperparam.H      = H;

% parameters
param.ksi    = ksi;
param.chi    = chi;
param.lambda = lambda;
param.kappa  = kappa;

% define data file name
file_name = [case_name,'.dat'];

% open data file
fileID = fopen(file_name,'w');

% define data format
formatSpec = '%+0.4f %0.4f %1.4f %.4f %.4f  \n';

% print global optimum in the data file
fprintf(fileID,'\n S_max   P_max  K01    mu1    mu2 \n');

for np1 = 1:Np1
    
	% update design variable
	f = X1(np1);
    
    for np2 = 1:Np2
        
        ngrid = np2 + (np1-1)*Np2;

        % print loop indicador
        if mod(ngrid,sqrt_Ngrid) == 0
            disp(['ngrid = ',num2str(ngrid)]);
        end
       
        % update design variable
        Omega = X2(np2);
        
        % define physical paramters vector
        param.f     = f;
        param.Omega = Omega;
                          
        % penalized objective function S(x)
        [S(ngrid),power(ngrid),K01(ngrid)] = ...
                    PiezoMagBeam_ObjFuncNoise(param,hyperparam);
        
        % update global maximum
        if S(ngrid) > S_opt
            
              S_opt =     S(ngrid);
          power_opt = power(ngrid);
            K01_opt =   K01(ngrid);
              X_opt = [f; Omega];
        end
        
        % print local values in data file
        fprintf(fileID,formatSpec,...
                S(ngrid),power(ngrid),K01(ngrid),X1(np1),X2(np2));
        
    end
end

% print global optimum in the data file
%fprintf(fileID,'\n');
fprintf(fileID,'\n\nGlobal maximum (ObjFunc power K01 f Omega):\n');
fprintf(fileID,formatSpec,S_opt,power_opt,K01_opt,X_opt(1),X_opt(2));

% close data file
fclose(fileID);

% display global optimum on the screen
fprintf('\n\nGlobal maximum (ObjFunc power K01 f Omega):\n');
fprintf(formatSpec,S_opt,power_opt,K01_opt,X_opt(1),X_opt(2));
fprintf('\n');

disp(' ');
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


% post processing
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');

% map to a 2d data structure
% (x = columns and y = lines)
    S = reshape(S    ,Np2,Np1);
power = reshape(power,Np2,Np1);
  K01 = reshape(K01  ,Np2,Np1);


% plot performance function countourmap
% ......................................................
gtitle = ' objective function contour map';
xlab   = ' excitation amplitude';
ylab   = ' excitation frequency';
xmin   = X_min(1);
xmax   = X_max(1);
ymin   = X_min(2);
ymax   = X_max(2);
gname  = [case_name,'__objfunc'];
flag   = 'eps';
fig1   = graph_contour_pnt(X1,X2,S'/S_opt,...
                           X_opt(1),X_opt(2),...
                           gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1);
% ......................................................


% plot mean power countourmap
% ......................................................
gtitle = ' mean power contour map';
xlab   = ' excitation amplitude';
ylab   = ' excitation frequency';
xmin   = X_min(1);
xmax   = X_max(1);
ymin   = X_min(2);
ymax   = X_max(2);
gname  = [case_name,'__mean_power'];
flag   = 'eps';
fig2   = graph_contour_pnt(X1,X2,power'/power_opt,...
                           X_opt(1),X_opt(2),...
                           gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

%close(fig2);
% ......................................................


% plot 0-1 identification parameter
% ......................................................
xlab   = 'excitation amplitude';
ylab   = 'excitation frequency';
label0 = 'regular';
label1 = 'chaos';
gtitle = '0-1 test for chaos';
xmin   = X_min(1);
xmax   = X_max(1);
ymin   = X_min(2);
ymax   = X_max(2);
zmin   = 0;
zmax   = 1;
gname  = [case_name,'__01test'];
flag   = 'eps';
fig3   = graph_binarymap(X1,X2,K01,...                          
                         X_opt(1),X_opt(2),...
                         gtitle,xlab,ylab,label0,label1,...
                         xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig3);
% ......................................................

toc
% -----------------------------------------------------------

% -----------------------------------------------------------------
%  HarvesterOpt_DirectSearch4var.m
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
%  The design variables are: (ksi,chi,lambda,kappa)
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
disp(' Piezo-Magneto-Elastic Harvester Optimizer          ')
disp(' (direct search)                                    ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Cunha (UERJ)                               ')
disp(' americo.cunha@uerj.br                              ')
disp(' ---------------------------------------------------')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'HarvesterOpt_DS_4var';

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

%ksi    = 0.01;  % mechanical damping ratio
%chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
%lambda = 0.05;  % dimensionless time constant reciprocal
%kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
f      = 0.083; % dimensionless excitation amplitude
Omega  = 0.8;   % dimensionless excitation frequency

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
Ndv = 4;

% number of points for design variables discretization
Np1 = 16;
Np2 = 16;
Np3 = 16;
Np4 = 16;

% number of points in the numerical grid
Ngrid = Np1*Np2*Np3*Np4;

% square root of Ngrid
sqrt_Ngrid = round(sqrt(Ngrid));

% penalty parameter
H = 10.0;

% limits for the design variables: (ksi,chi,lambda,kappa)
X_min = [0.01; 0.050; 0.050; 0.50];
X_max = [0.05; 0.200; 0.200; 1.50];

% numerical grid for domain discretization
X1 = linspace(X_min(1),X_max(1),Np1);
X2 = linspace(X_min(2),X_max(2),Np2);
X3 = linspace(X_min(3),X_max(3),Np3);
X4 = linspace(X_min(4),X_max(4),Np4);

% preallocate memory for PerfFunc
S = zeros(Ngrid,1);

% preallocate memory for power function
power = zeros(Ngrid,1);

% preallocate memory for 0-1 classifier
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
param.f     = f;
param.Omega = Omega;

% define data file name
file_name = [case_name,'.dat'];

% open data file
fileID = fopen(file_name,'w');

% define data format
formatSpec = '%+0.4f %0.4f %1.4f %.4f %.4f %.4f %.4f \n';

% print global optimum in the data file
fprintf(fileID,'\n S_max   P_max  K01    mu1    mu2    mu3    mu4 \n');

for np1 = 1:Np1
    
	% update design variable
    ksi = X1(np1);
    
    for np2 = 1:Np2
        
        % update design variable
        chi = X2(np2);
            
        for np3 = 1:Np3

        % update design variable
        lambda = X3(np3);
            
            for np4 = 1:Np4
                
                ngrid = np4 + (np3-1)*Np4 + ...
                              (np2-1)*Np4*Np3 + ...
                              (np1-1)*Np4*Np3*Np2;

                % print loop indicador
                if mod(ngrid,sqrt_Ngrid) == 0
                    disp(['ngrid = ',num2str(ngrid)]);
                end

                % update design variable
                kappa  = X4(np4);

                % define physical paramters vector
                param.ksi    = ksi;
                param.chi    = chi;
                param.lambda = lambda;
                param.kappa  = kappa;

                % penalized objective function S(x)
                [S(ngrid),power(ngrid),K01(ngrid)] = ...
                            PiezoMagBeam_ObjFunc(param,hyperparam);

                % update global maximum
                if S(ngrid) > S_opt

                      S_opt =     S(ngrid);
                  power_opt = power(ngrid);
                    K01_opt =   K01(ngrid);
                      X_opt = [ksi; chi; lambda; kappa];
                end
                
                % print local values in data file
                fprintf(fileID,formatSpec,...
                        S(ngrid),power(ngrid),K01(ngrid),...
                        X1(np1),X2(np2),X3(np3),X4(np4));
            end
        end
    end
end

% print global optimum in the data file
fprintf(fileID,'\n\nGlobal maximum (ObjFunc power K01 ksi chi lambda kappa):\n');
fprintf(fileID,formatSpec,...
        S_opt,power_opt,K01_opt,X_opt(1),X_opt(2),X_opt(3),X_opt(4));

% close data file
fclose(fileID);

% display global optimum on the screen
fprintf('\n\nGlobal maximum (ObjFunc power K01 ksi chi lambda kappa):\n');
fprintf(formatSpec,...
        S_opt,power_opt,K01_opt,X_opt(1),X_opt(2),X_opt(3),X_opt(4));
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

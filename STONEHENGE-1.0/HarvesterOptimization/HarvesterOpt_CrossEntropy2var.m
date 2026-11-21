% -----------------------------------------------------------------
%  HarvesterOpt_CrossEntropy2var.m
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
%  The solution strategy is based on the cross-entropy method,
%  a stochastic metaheuristic that 'transforms' the optimization 
%  problem into a rare event estimation problem.
%  
%  The design variables are: (f,omega)
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
disp(' Piezo-Magneto-Elastic Harvester Optimization       ')
disp(' (cross-entropy method)                             ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Cunha (UERJ)                               ')
disp(' americo.cunha@uerj.br                              ')
disp(' ---------------------------------------------------')
% -----------------------------------------------------------



% simulation information
% -----------------------------------------------------------
case_name = 'HarvesterOpt_CE_2var';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% define physical parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining physical parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');

ksi    = 0.01;  % mechanical damping ratio
chi    = 0.05;  % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;  % dimensionless time constant reciprocal
kappa  = 0.5;   % dimensionless piezoeletric coupling term (eletrical)
%f      = 0.115; % dimensionless excitation amplitude
%Omega  = 0.8;   % dimensionless excitation frequency

x0    = 1.0; % dimensionless initial displacement
xdot0 = 0.0;  % dimensionless initial velocity
v0    = 0.0;  % dimensionless initial voltage

toc
% -----------------------------------------------------------


% 0-1 test for chaos parameters
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- defining model parameters --- ');
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


% solve optimization problem via cross-entropy method
%---------------------------------------------------------------
tic
disp(' '); 
disp(' --- optimization via cross-entropy method --- ');
disp(' ');
disp('    ... ');
disp(' ');

% number of design variables
Ndv = 2;

% maximum number of iterations
maxiter = 100;

% cross-entropy tolerance
tolCE = 1.0e-3;

% smoothing parameter (0 < alpha <= 1) 
% -- set alpha = 1 for no smoothing --
alpha = 0.7;
one_minus_alpha = 1.0 - alpha;

% dynamic smoothing parameters
% (0.8 <= beta <= 0.99)
% (q is a interger between 5 and 10)
beta = 0.8;
q = 5;

% penalty parameter
H = 10.0;

% elite samples percentage
rho = 0.1;

% number of cross-entropy samples
NCE = 50;

% elite samples index
Nelite = NCE - ceil(rho*NCE) + 1;

% initialize level counter
t = 0;

% preallocate memory for PerfFunc
S = zeros(NCE,maxiter);

% preallocate memory for power function
power = zeros(NCE,maxiter);

% preallocate memory for 0-1 classifier
K01 = zeros(NCE,maxiter);

% preallocate memory for design variables samples
%X = zeros(NCE*maxiter,Ndv);
X1 = zeros(NCE,maxiter);
X2 = zeros(NCE,maxiter);

% initialize ObjFunc maximum
S_opt = -Inf;

% initialize power maximum
power_opt = -Inf;

% initialize K01 maximum
K01_opt = -Inf;

% initialize maximum point
X_opt = zeros(Ndv,1);

% limits for the design variables: (f,Omega)
X_min = [0.08; 0.75];
X_max = [0.10; 0.85];

% initialize mean and std dev vectors
   mu = X_min + (X_max-X_min).*rand(Ndv,1);
sigma =       5*(X_max-X_min);

% initialize old mean and std dev vectors
   mu0 = mu;
sigma0 = sigma;

% inf-norm of sigma
sigma_max = norm(sigma,Inf);

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
formatSpec1 = '%03d %+0.4f %0.4f %1.4f %.4f %.4f %.4f %.4f \n';
formatSpec2 = '%+0.4f %0.4f %1.4f %.4f %.4f \n';

% print global optimum in the data file
fprintf(fileID,'\n t   S_max   P_max  K01    mu1    mu2   sigma1 sigma2\n');

% display extreme values on screen
fprintf('\n t   S_max   P_max  K01    mu1    mu2   sigma1 sigma2\n');

while sigma_max > tolCE && t < maxiter
        
        % update level counter
        t = t + 1;
        
        % limit vectors for truncated Gaussian
        supp_p_l = ((X_min - mu)./sigma)*ones(1,NCE);
        supp_p_h = ((X_max - mu)./sigma)*ones(1,NCE);
        %supp_X1_l = ((X_min(1) - mu(1))./sigma(1))*ones(NCE,1);
        %supp_X1_h = ((X_max(1) - mu(1))./sigma(1))*ones(NCE,1);
        %supp_X2_l = ((X_min(2) - mu(2))./sigma(2))*ones(NCE,1);
        %supp_X2_h = ((X_max(2) - mu(2))./sigma(2))*ones(NCE,1);
        
        % generate samples from truncated normal distribution
        X1(:,t) = mu(1) + trandn(supp_p_l(1,:),supp_p_h(1,:))*sigma(1);
        X2(:,t) = mu(2) + trandn(supp_p_l(2,:),supp_p_h(2,:))*sigma(2);
        %X1(:,t) = mu(1) + trandn(supp_X1_l,supp_X1_h)*sigma(1);
        %X2(:,t) = mu(2) + trandn(supp_X2_l,supp_X2_h)*sigma(2);
        
        % evaluate objective function at the samples
        for n=1:NCE
        
            % update parameters
            param.f     = X1(n,t);
            param.Omega = X2(n,t);
                          
            % penalized objective function S(x)
            [S(n,t),power(n,t),K01(n,t)] = ...
                    PiezoMagBeam_ObjFunc(param,hyperparam);
        end
        
        % sort objective function evaluations
        [S_sort,I_sort] = sort(S(:,t));
        
        % update mean value
        mu(1) = mean(X1(I_sort(Nelite:end),t));
        mu(2) = mean(X2(I_sort(Nelite:end),t));

        % update standard deviation
        sigma(1) = std(X1(I_sort(Nelite:end),t));
        sigma(2) = std(X2(I_sort(Nelite:end),t));
        
        % estimator for objective function maximum
          S_max = S_sort(Nelite);
          P_max = power(I_sort(Nelite),t);
        K01_max =   K01(I_sort(Nelite),t);
        
        % smoothing
        mu = alpha*mu + one_minus_alpha*mu0;
        
        % update dynamics smoothing parameter
        beta_t = beta*(1 - (1-1/t)^q);

        % dynamic smoothing
        sigma = beta_t*sigma + (1-beta_t)*sigma0;
        
        % save old mean
        mu0 = mu;

        % save old std dev
        sigma0 = sigma;
        
        % inf-norm of sigma
        sigma_max = norm(sigma,Inf);
        
        % update global maximum
        if S_max > S_opt
             S_opt = S_max;
         power_opt = P_max;
           K01_opt = K01_max;
          X_opt(1) = X1(I_sort(Nelite),t);
          X_opt(2) = X2(I_sort(Nelite),t);
        end
            
        % print local extreme values in a data file
        fprintf(fileID,formatSpec1,...
                t,S_max,P_max,K01_max,mu(1),mu(2),sigma(1),sigma(2));
        
        % print local extreme values on screen
        fprintf(formatSpec1,...
                t,S_max,P_max,K01_max,mu(1),mu(2),sigma(1),sigma(2));
end

% print global optimum in the data file
fprintf(fileID,'\n\nGlobal maximum (ObjFunc power K01 f Omega):\n');
fprintf(fileID,formatSpec2,...
        S_opt,power_opt,K01_opt,X_opt(1),X_opt(2));
fprintf('\n');

% close data file
fclose(fileID);


% display global optimum on the screen
fprintf('\n\nGlobal maximum (ObjFunc power K01 f Omega):\n');
fprintf(formatSpec2,...
        S_opt,power_opt,K01_opt,X_opt(1),X_opt(2));


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




% animate cross-entropy iteration
% ...........................................................
vtitle = 'cross-entropy method';
xlab   = 'excitation amplitude';
ylab   = 'excitation frequency';
xmin   = X_min(1);
xmax   = X_max(1);
ymin   = X_min(2);
ymax   = X_max(2);
xref = 0.0994; % obtained via direct search with a 256 x 256 grid
yref = 0.7771; % obtained via direct search with a 256 x 256 grid
vname  = [num2str(case_name),'__animation'];

plot_ce_animation(X1,X2,(1:t),xref,yref,...
                   vtitle,xlab,ylab,xmin,xmax,ymin,ymax);
% ...........................................................




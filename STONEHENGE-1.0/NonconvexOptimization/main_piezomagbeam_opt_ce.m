
% -----------------------------------------------------------------
%  main_piezomagbeam_opt_ce.m
%
%  This is the main file for a program which solve an optimization
%  problem, via cross-entropy method, that seaks to maximize 
%  output power of a piezo-magneto-elastic beam, which evolves 
%  according to
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
%  A. Cunha Jr
%  Enhancing the performance of a bistable energy harvesting
%  device via the cross-entropy method
%  2020
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: March 6, 2017
% -----------------------------------------------------------------


clc
clear all
close all


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' Piezo-Magneto-Elastic Beam Dynamics                ')
disp(' (optimization via cross-entropy)                   ')
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
f      = 0.083; % dimensionless excitation amplitude
Omega  = 0.08;  % dimensionless excitation frequency

x0    = 1.0; % dimensionless initial displacement
xdot0 = 0.0;  % dimensionless initial velocity
v0    = 0.0;  % dimensionless initial voltage
% -----------------------------------------------------------


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


% ODE solver parameters
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
% -----------------------------------------------------------


% solve optimization problem via cross-entropy method
%---------------------------------------------------------------
tic

disp(' '); 
disp(' --- optimization problem via cross-entropy --- ');
disp(' ');
disp('    ... ');
disp(' ');


% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setDefaultStream(rng_stream); % Matlab 2009
%RandStream.setGlobalStream(rng_stream); % Matlab 2013

% maximum number of iterations
maxiter = 100;

% cross-entropy tolerance
tolCE = 1.0e-3;

% fixed smoothing parameter
% (0 <= alpha <= 1) typically 0.7
% alpha = 1 for no smoothing
alpha = 0.7;
one_minus_alpha = 1.0 - alpha;

% dynamic smoothing parameters
% (0.8 <= beta <= 0.99)
% (q is a interger between 5 and 10)
beta = 0.8;
q = 5;

% penalty parameter
Hpenalty = 10.0;

% elite samples percentage
rho = 0.1;

% number of cross-entropy samples
NCE = 25;

% elite samples index
Nelite = NCE - ceil(rho*NCE) + 1;
%Nelite = NCE -10+1;

% initialize level counter
t = 0;

% preallocate memory for PerfFunc
S = zeros(NCE,maxiter);

% preallocate memory for power function
power = zeros(NCE,maxiter);

% preallocate memory for 0-1 classifier
K01 = zeros(NCE,maxiter);

% preallocate memory for p1 samples (excitation amplitude)
p1 = zeros(NCE,maxiter);

% preallocate memory for p2 samples (excitation frequency)
p2 = zeros(NCE,maxiter);

% initialize PerfFunc maximum value
S_max_overall = -Inf;

% initialize power maximum value
power_max_overall = -Inf;

% initialize K01 maximum value
K01_max_overall = -Inf;

% initialize PerfFunc maximum point
p1_max_overall = 0.0;
p2_max_overall = 0.0;

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

% support for initial parameter vector
pmin = [p1_min p2_min];
pmax = [p1_max p2_max];

% initialize mean and std dev vectors
   mu = pmin + (pmax-pmin).*rand(1,2);
sigma =      5*(pmax-pmin);

% initialize old mean and std dev vectors
   mu0 = mu;
sigma0 = sigma;

% simulation information
%case_name  = ['piezomagbeam_opt_ce_',...
%                'penal_',num2str(Hpenalty),'_',...
%                  'NCE_',num2str(NCE),'_',...
%              'maxiter_',num2str(maxiter)];

% simulation information
case_name  = ['piezomagbeam_optdesign_ce_',...
                'penal_',num2str(Hpenalty),'_',...
                  'NCE_',num2str(NCE),'_',...
              'maxiter_',num2str(maxiter)];

% define data file name
file_name = [case_name,'.dat'];

% open data file
fileID = fopen(file_name,'w');

% define data format
formatSpec1 = '%03d %+0.4E %+0.4E %+0.4E %0.4E %0.4E \n';

% display extreme values on screen
%fprintf('\n t    S_max       mu1         mu2         sigma1      sigma2\n');

while norm(sigma,Inf) > tolCE && t < maxiter
        
        % update level counter
        t = t+1;
        
        % generate samples from normal distribution
        %p1(:,t) = mu(1,1) + randn(NCE,1)*sigma(1,1);
        %p2(:,t) = mu(1,2) + randn(NCE,1)*sigma(1,2);
        
        % limit vectors for truncated Gaussian
        supp_p1_l = ((pmin(1,1) - mu(1,1))./sigma(1,1))*ones(NCE,1);
        supp_p1_h = ((pmax(1,1) - mu(1,1))./sigma(1,1))*ones(NCE,1);
        supp_p2_l = ((pmin(1,2) - mu(1,2))./sigma(1,2))*ones(NCE,1);
        supp_p2_h = ((pmax(1,2) - mu(1,2))./sigma(1,2))*ones(NCE,1);
        
        % generate samples from truncated normal distribution
        p1(:,t) = mu(1,1) + trandn(supp_p1_l,supp_p1_h)*sigma(1,1);
        p2(:,t) = mu(1,2) + trandn(supp_p2_l,supp_p2_h)*sigma(1,2);
        
        % evaluate performance function at the samples
        for n=1:NCE
        
            % update parameters
            %f     = p1(n,t);
            %Omega = p2(n,t);
            chi   = p1(n,t);
            kappa = p2(n,t);
        
            % define physical paramters vector
            phys_param = [ksi chi f Omega lambda kappa];
                          
            % penalized performance function S(x)
            [S(n,t),power(n,t),K01(n,t)] = ...
              piezomagbeam_opt_PerfFunc(phys_param,tspan,IC,...
                                        cmin,cmax,Nc,tol01,OSflag,Hpenalty,...
                                        p1_min,p1_max,p2_min,p2_max);
        end
        
        % sort performance function evaluations
        [S_sort,I_sort] = sort(S(:,t));
        
        % update mean value
        mu(1,1) = mean(p1(I_sort(Nelite:end),t));
        mu(1,2) = mean(p2(I_sort(Nelite:end),t));

        % update standard deviation
        sigma(1,1) = std(p1(I_sort(Nelite:end),t));
        sigma(1,2) = std(p2(I_sort(Nelite:end),t));
        
        % estimator for PerfFunc maximum
        S_max = S_sort(Nelite);
        
        % fixed smoothing
           mu = alpha*mu    + one_minus_alpha*mu0;
        %sigma = alpha*sigma + one_minus_alpha*sigma0;
        
        % update dynamics smoothing parameter
        beta_t = beta*(1 - (1-1/t)^q);

        % dynamic smoothing
        sigma = beta_t*sigma + (1-beta_t)*sigma0;
        
        % save old mean
        mu0 = mu;

        % save old std dev
        sigma0 = sigma;
        
        % update global maximum
        if S_max > S_max_overall
             S_max_overall = S_max;
         power_max_overall = power(I_sort(Nelite),t);
           K01_max_overall =   K01(I_sort(Nelite),t);
            p1_max_overall =    p1(I_sort(Nelite),t);
            p2_max_overall =    p2(I_sort(Nelite),t);
        end
        
        % print local extreme values on screen
        fprintf(formatSpec1,...
                t,S_max,mu(1),mu(2),sigma(1),sigma(2));
            
        % print local extreme values in a data file
        fprintf(fileID,formatSpec1,...
                t,S_max,mu(1),mu(2),sigma(1),sigma(2));
end

% define data format
formatSpec2 = '%.4f %.4f %+0.4E %0.4E %1.4f \n';

% display global extreme values on screen
fprintf('\n\nGlobal maximum (f Omega PerfFunc power K01):\n\n');

fprintf(formatSpec2,...
        p1_max_overall,p2_max_overall,...
        S_max_overall,power_max_overall,K01_max_overall);

% display global extreme values in a data file
fprintf(fileID,'\n\nGlobal maximum (f Omega PerfFunc power K01):\n\n');
fprintf(fileID,formatSpec2,...
        p1_max_overall,p2_max_overall,...
        S_max_overall,power_max_overall,K01_max_overall);

fprintf('\n');

% close data file
fclose(fileID);

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




% animate cross-entropy iteration
% ...........................................................
%global max: 0.0994 0.7771 +3.4378E-01 3.4378E-01 0.0393 
%global min: 0.0802 0.7771 -8.9856E+00 6.8530E-03 0.9992 
vtitle = 'cross-entropy method';
%xlab   = 'excitation amplitude';
%ylab   = 'excitation frequency';
xlab   = ' mechanical coupling';
ylab   = ' electrical coupling';
xmin   = p1_min;
xmax   = p1_max;
ymin   = p2_min;
ymax   = p2_max;
xref = 0.0994;
yref = 0.7771;
vname  = [num2str(case_name),'__animation'];

plot_ce_animation(p1,p2,(1:t),xref,yref,...
                        vtitle,xlab,ylab,xmin,xmax,ymin,ymax);
%mov = plot_ce_animation_video(p1,p2,(1:t),...
%                              xref,yref,vtitle,vname,...
%                              xlab,ylab,xmin,xmax,ymin,ymax);
% ...........................................................




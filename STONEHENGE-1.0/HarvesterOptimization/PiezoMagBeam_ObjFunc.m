% -----------------------------------------------------------------
%  PiezoMagBeam_ObjFunc.m
% ----------------------------------------------------------------- 
%  Objective function for the optimization problem.
% ----------------------------------------------------------------- 
%  programmer: Americo Cunha
%              americo.cunha@uerj.br
%
%  last update: October 23, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [S,power,K01] = PiezoMagBeam_ObjFunc(param,hyperparam)

% hyper parameters
tspan  = hyperparam.tspan;
IC     = hyperparam.IC;
cmin   = hyperparam.cmin;
cmax   = hyperparam.cmax;
Nc     = hyperparam.Nc;
tol01  = hyperparam.tol01;
OSflag = hyperparam.OSflag;
H      = hyperparam.H;

% parameters
%f      = param.f;
%Omega  = param.Omega;
%ksi    = param.ksi;
%chi    = param.chi;
lambda = param.lambda;
%kappa  = param.kappa;


% integrate the dynamical system with RK-45
[time,y] = ode45(@(t,y)PiezoMagBeam_RHS(t,y,param),tspan,IC);

% number of time steps
Ndt = length(time);

% number of steps for steady state
Nss = round(0.5*Ndt);

% number of jumps
Njump = round((Ndt-Nss)/(0.05*Ndt));

% temporal interval of analysis
T = time(end)-time(Nss);

% compute the mean output power
power = (lambda/T)*trapz(time(Nss:end),y(Nss:end,3).^2);
        
% apply 0-1 test for chaos in voltage time series
K01 = test01chaos(y(Nss:Njump:end,3),cmin,cmax,Nc,OSflag);

% penalization function
Penalty = max(K01 - tol01,0);
                          
% penalized objective function
S = power - H*Penalty;

end
% -----------------------------------------------------------------

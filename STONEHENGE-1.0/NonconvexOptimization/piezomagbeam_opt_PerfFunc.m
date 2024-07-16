
% -----------------------------------------------------------------
%  piezomagbeam_opt_PerfFunc.m
%
%  This function computes the performance function for
%  cross entropy method.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: March 8, 2017
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [S,power,K01] = piezomagbeam_opt_PerfFunc(phys_param,...
                                                   tspan,IC,...
                                                   cmin,cmax,Nc,tol01,...
                                                   OSflag,Hpenalty,...
                                                   p1_min,p1_max,...
                                                   p2_min,p2_max)

% physical parameters
%ksi    = phys_param(1);
%chi    = phys_param(2);
f      = phys_param(3);
Omega  = phys_param(4);
%lambda = phys_param(5);
%kappa  = phys_param(6);
%x0     = phys_param(7);
%xdot0  = phys_param(8);
%v0     = phys_param(9);

% integrate the dynamical system with RK-45
[time,y] = ode45(@(t,y)piezomagbeam(t,y,phys_param),tspan,IC);

% number of time steps
Ndt = length(time);

% number of steps for steady state
Nss = round(0.5*Ndt);

% number of jumps
Njump = round((Ndt-Nss)/(0.05*Ndt));

% temporal interval of analysis
T = time(end)-time(Nss);

% compute the output power
power = (1/T)*trapz(time(Nss:end),y(Nss:end,3).^2);
        
% apply 0-1 test for chaos in voltage time series
K01 = test01chaos(y(Nss:Njump:end,3),cmin,cmax,Nc,OSflag);
                          
% penalized performance function
S = power - ...
    Hpenalty*max((   K01 -  tol01),0) + ...
    Hpenalty*max((     f - p1_max),0) + ...
    Hpenalty*max((    -f + p1_min),0) + ...
    Hpenalty*max(( Omega - p2_max),0) + ...
    Hpenalty*max((-Omega + p2_min),0);

end
% -----------------------------------------------------------------

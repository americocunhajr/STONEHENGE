% -----------------------------------------------------------------
%  piezomagbeam_asymmetric.m
%
%  This function defines the right hand side of the follwoing
%  nonlinear system of ordinary differential equations
%
%    dy1/dt = y2
%    dy2/dt = -2*ksi.y2 + 0.5*y1(1+2.delta.y1-y1^2) + (1+beta.|y1|).chi.y3 + f.cos(Omega.t)
%    dy3/dt = - lambda.y3 - (1+beta.|y1|).kappa.y2,
%
%  that models the dynamics of a piezo-magneto-elastic beam.
%  
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Jul 05, 2021
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [time,Y] = piezomagbeam_asymmetric(X,IC,tspan)
% check number of arguments
     if nargin > 3
         error('Too many inputs.')
     elseif nargin < 3
         error('Too few inputs.')
     end 
     
% check number of parameters 
    if numel(fieldnames(X)) > 9 
        error('Too many parameters inputs.')
    elseif numel(fieldnames(X)) < 9
         error('Too few parameters inputs.')
    end 
     
% physical parameters
    ksi    = X.ksi;
    chi    = X.chi;
    lambda = X.lambda;
    kappa  = X.kappa;
    f      = X.f;
    Omega  = X.Omega; 
    beta   = X.beta;
    delta  = X.delta;
    phi    = X.phi;
    p      = 0.59;
    
% ODE solver optional parameters
    opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
    
% ODE right hand side (Z = [z; zdot])
    dYdt = @(t,y)[y(2);
                  -2*ksi*y(2) + 0.5*y(1)*(1.0+2*delta*y(1)-y(1)^2) + (1+beta*abs(y(1)))*chi*y(3) + p*sin(phi*pi/180) + f*cos(Omega*t);
                  -lambda*y(3) - (1+beta*abs(y(1)))*kappa*y(2)];
              
% ODE solver Runge-Kutta45
    [time,Y] = ode45(dYdt,tspan,IC,opt);

end
% -----------------------------------------------------------------
%  HarvesterColoredNoise.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% -----------------------------------------------------------------
%  This function defines the right hand side of the following
%  system of nonlinear ordinary differential equations
%
%  dy1dt = y2
%  dy2dt = -2 ksi y2 + 0.5y1(1-y1^2) + chi y3 + fcos(Omega t) + noise
%  dy3dt = - lambda y3 - kappa y2,
%
%  that models the nonlinear stochastic dynamics of a bistable 
%  piezo-magneto-elastic energy harvester.
%  
%  Input:
%  t     - time (s)
%  y     - state space vector
%  param - physical parameters structure
%  
%  Output:
%  dydt - right hand side function
% ----------------------------------------------------------------- 

% -----------------------------------------------------------------
function dydt = HarvesterColoredNoise(t,y,param)



% physical parameters
ksi    = param.ksi;
chi    = param.chi;
lambda = param.lambda;
kappa  = param.kappa;
f      = param.f;
Omega  = param.Omega;

% noise time-series
noise = param.noise;

% temporal mesh with the instants where the noise is computed
time_noise = param.time_noise;

% evalute the noise at the time instant t
noise_t = interp1(time_noise,noise,t,'linear');


% y = [x dxdt v] is the state vector
x    = y(1);
dxdt = y(2);
v    = y(3);

% dydt = [dxdt d2xdt2 dvdt] is the evolution law
d2xdt2 = -2*ksi*dxdt + 0.5*x*(1-x^2) + chi*v + f*cos(Omega*t) + noise_t;
dvdt   = -lambda*v - kappa*dxdt;
dydt   = [dxdt; d2xdt2; dvdt];

end
% -----------------------------------------------------------------

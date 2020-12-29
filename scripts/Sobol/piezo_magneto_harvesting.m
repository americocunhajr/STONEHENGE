% -----------------------------------------------------------------
%  piezo_magneto_harvesting.m
%
%  This function computes the power of a piezo-magneto-elastic
%  beam.
% ----------------------------------------------------------------- 
%  programmer: João Pedro C V Norenberg
%              jpcvalese@gmail.com
%
%  last update: Dez 13, 2019
% -----------------------------------------------------------------
%  Reference:
%  
%  A. Erturk, J. Hoffmann, and D. J. Inman
%  A piezomagnetoelastic structure for broadband vibration
%  energy harvesting
%  Applied Physics Letters
%  vol. 94 pp. 254102, 2009

function mean_power = piezo_magneto_harvesting(X)

% model parameters
    ksi     = X(1);
    chi     = X(2);
    lambda  = X(3);
    kappa   = X(4);
    f       = X(5);
    Omega   = X(6);

% initial time of analysis (s)
    t0 = 0.0;

% initial time of analysis (s)
    t1 = 2000.0;

% increment time
    Ndt = 0.01;

% time interval of analysis
    tspan = t0:Ndt:t1;

% initial displacement
    y0 = 1.0;

% initial velocity
    ydot0 = 0.0;

% initial voltage
    v0 = 0.0;
    
 % initial conditions
IC = [y0 ydot0 v0]';

% ODE solver optional parameters
opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

% ODE right hand side (Z = [z; zdot])
dYdt = @(t,y) [y(2);
               -2.*ksi.*y(2) + 0.5.*y(1).*(1.0-y(1).^2) + chi.*y(3) + f.*cos(Omega.*t);
               -lambda.*y(3) - kappa.*y(2)];
               

% ODE solver Runge-Kutta45
[time,Y] = ode45(dYdt,tspan,IC,opt);

begin1 = round(time(end)*2/3);
% temporal voltage
Qvolt = Y(begin1:end,3);
% temporal interval of analysis
T = time(end)-time(begin1);

% output power
power = lambda.*Qvolt.^2;
mean_power = (1/T).*trapz(time(begin1:end),power);

end
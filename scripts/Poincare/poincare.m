% -----------------------------------------------------------------
%  poincare.m
% ----------------------------------------------------------------- 
%  This function compute the Poincare section.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%
%  last update: Dec 20, 2020
% -----------------------------------------------------------------

function [poincare_disp, poincare_velo] = poincare(f,beta)
    % check number of arguments
    if nargin > 2
        error('Too many inputs.')
    elseif nargin < 2
        error('Too few inputs.')
    end 
    dt = 0.0005;

    ksi    = 0.01;          % mechanical damping ratio
    chi    = 0.05;          % dimensionless piezoeletric coupling term (mechanical)
    lambda = 0.05;          % dimensionless time constant reciprocal
    kappa  = 0.5;           % dimensionless piezoeletric coupling term (eletrical)
    Omega  = 0.8;           % dimensionless excitation frequency
    IC     = [1,0,0];       % initial conditions 
    
    % defining equation of motion
    func = @(t,y) [y(2);
        -2.*ksi.*y(2) + 0.5.*y(1).*(1.0-y(1).^2) + (1+beta*abs(y(1)))*chi.*y(3) + f.*cos(Omega.*t);
        -lambda.*y(3) - (1+beta*abs(y(1)))*kappa.*y(2)];

    % period of a forcing cycle
    T = 2*pi/Omega;

    % number of forcing cycles
    Nf =  5000;

    % initial dimensionless time
    t0 = 0.0;

    % final dimensionless time
    t1 = t0 + Nf*T;

    % number of samples per forcing cycle
    Nsamp = 10000;

    % time series sampling points
    tspan = t0:(T/Nsamp):t1;

    % ODE solver optional parameters
    opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

    % ODE solver Runge-Kutta45
    [time,y] = ode45(func,tspan,IC,opt);

    % time series of dimensionless displacement
    Qdisp = y(:,1);

    % time series of dimensionless velocity
    Qvelo = y(:,2);

    % time series of dimensionless voltage
    Qvolt = y(:,3);

    % number of dimensionless time steps
    Ndt = length(time);

    % number of steps for steady state
    Nss = round(0.99*Ndt);

    % number of steps to initiates Poincare map
    Npm = round(0.10*Ndt);

    % Poincare maps
    poincare_disp = Qdisp(Npm:Nsamp:Ndt);
    poincare_velo = Qvelo(Npm:Nsamp:Ndt);
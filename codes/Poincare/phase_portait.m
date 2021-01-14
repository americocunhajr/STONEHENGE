% -----------------------------------------------------------------
%  phase_portait.m
% ----------------------------------------------------------------- 
%  This function compute the phase portait.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Dec 20, 2020
% -----------------------------------------------------------------

function [disp,vel] = phase_portait(f,beta)

    % check number of arguments
    if nargin > 2
        error('Too many inputs.')
    elseif nargin < 2
        error('Too few inputs.')
    end 
    
    % input parameters
    ksi       = 0.01;           % damping ratio
    chi       = 0.05;           % piezoelectric coupling (mechanical)
    lambda    = 0.05;           % reciprocal time constant
    kappa     = 0.5;            % piezoelectric coupling (electrical)
    Omega     = 0.8;            % frequency of external force
    x0        = [1,0,0];        % initial condition

    % number of forcing cycles
    Nf =  10000;

    % initial dimensionless time
    t0 = 0.0;

    t1 = t0 + Nf;

    tspan = t0:0.001:t1;
    
    func = @(t,y) [y(2);
        -2.*ksi.*y(2) + 0.5.*y(1).*(1.0-y(1).^2) + (1+beta*abs(y(1)))*chi.*y(3) + f.*cos(Omega.*t);
        -lambda.*y(3) - (1+beta*abs(y(1)))*kappa.*y(2)];

    opt1 = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

    % ODE solver Runge-Kutta45
    [~,Y1] = ode45(func,tspan,x0,opt1);
    disp = Y1(round(0.5*end):end,1);
    vel = Y1(round(0.5*end):end,2);
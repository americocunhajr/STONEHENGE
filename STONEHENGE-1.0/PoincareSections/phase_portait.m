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

function [disp,vel] = phase_portait(Xpar,IC)

    % check number of arguments
    if nargin > 2
        error('Too many inputs.')
    elseif nargin < 2
        error('Too few inputs.')
    end 

    % number of forcing cycles
    Nf =  10000;

    % initial dimensionless time
    t0 = 0.0;

    t1 = t0 + Nf;

    tspan = t0:0.001:t1;
    
    % ODE solver Runge-Kutta45
    [~,Y1] = piezomagbeam_asymmetric(Xpar,IC,tspan);

    disp = Y1(round(0.5*end):end,1);
    vel = Y1(round(0.5*end):end,2);
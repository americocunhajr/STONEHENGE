% -----------------------------------------------------------------
%  poincare.m
% ----------------------------------------------------------------- 
%  This function compute the Poincare section.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Dec 20, 2020
% -----------------------------------------------------------------

function [poincare_disp, poincare_velo] = poincare(Xpar,IC)
    % check number of arguments
    if nargin > 2
        error('Too many inputs.')
    elseif nargin < 2
        error('Too few inputs.')
    end 
    
    Omega  = Xpar.Omega;           % dimensionless excitation frequency

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

    % ODE solver Runge-Kutta45
    [time,y] = piezomagbeam_asymmetric(Xpar,IC,tspan);

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
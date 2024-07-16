% -----------------------------------------------------------------
%  bifurcation_3d.m
% ----------------------------------------------------------------- 
%  This function defines the bifurcation diagram for amplitude of
%  excitation as control parameter.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Oct 20, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------

function [bifurc_inf] = bifurcation_3d(var_input,name_file)

    % check number of arguments
    if nargin > 2
        error('Too many inputs.')
    elseif nargin < 2
        error('Too few inputs.')
    end 
    
    % check number of field in struct
    if length(struct2cell(var_input)) < 5 
        error('Too few inputs in struct')
    elseif length(struct2cell(var_input)) > 5
        error('Too many inputs in struct')
    end
    
    % check number of parameters in field struct X_params
    if length(struct2cell(var_input.X_params)) < 7 
        error('Too few inputs in parameters model')
    elseif length(struct2cell(var_input.X_params)) > 7
        error('Too many inputs in parameters model')
    end
    
    % physical parameters                
    ksi    = var_input.X_params.ksi;    % mechanical damping ratio
    chi    = var_input.X_params.chi;    % piezoeletric coupling term (mechanical)
    lambda = var_input.X_params.lambda; % time constant reciprocal
    kappa  = var_input.X_params.kappa;  % piezoeletric coupling term (eletrical)
    beta   = var_input.X_params.beta;   % nonlinear electromechanical coupling term
    delta  = var_input.X_params.delta;  % asymmetric coefficient of potential energy
    phi    = var_input.X_params.phi;    % bias angle
    

    % bifurcation control parameters
    Omega_int = var_input.Par1_rang.Omega_rang; % Omega interval
    num_omega = var_input.N1_rang.N_omega;      % Size of Omega interval
    f_int     = var_input.Par2_rang.f_int;      % f interval
    int_param = var_input.N2_rang.int_param;    % Increment of f interval

    % range of frequency of excitation
    Omega_rang = linspace(min(Omega_int),max(Omega_int),num_omega);

    % bifurcation sweeping up f
    disp(' ')
    disp(' ... ')
    disp(' ---- bifurcation processing start --- ')
    disp(' ... ')

    % loop for Omega bifurcation
    for k = 1:length(Omega_rang)

    % Omega value in step 
    Omega = Omega_rang(k);

    % vector of sweeping up f
    f_rang_up = min(f_int):int_param:max(f_int);

    % initial condition
    x0 = [1,0,0];

    % loop for bifurcation processing 1 (sweeping up)
    for i = 1:length(f_rang_up)

        % f value in step
        f = f_rang_up(i);

        % struc variable physical parameters
        X_var.ksi    = ksi;
        X_var.chi    = chi;
        X_var.lambda = lambda;
        X_var.kappa  = kappa;
        X_var.beta   = beta;
        X_var.f      = f;
        X_var.Omega  = Omega;
        X_var.delta  = delta;
        X_var.phi    = phi;

        % period of a forcing cycle
        T = 2*pi/Omega;

        % number of forcing cycles
        Nf = 6000;

        % initial dimensionless time
        t0 = 0.0;

        % final dimensionless time
        t1 = t0 + Nf*T;

        % number of samples per forcing cycle
        Nsamp = 10;

        % time series sampling points
        tspan = t0:(T/Nsamp):t1;

        % cut transient regime
        Nss = round(length(tspan)*0.7);
        Npp = length(tspan);

        % compute solve differential equation
        [~,Y1] = piezomagbeam_asymmetric(X_var,x0,tspan);

        % vector of bifurcation point
        Bifurc_up(:,i,k) = Y1(Nss:Nsamp:Npp,3);

        % rewrite initial condition
        x0 = [Y1(end,1) Y1(end,2) Y1(end,3)];
    end

    % vector of sweeping down f
    f_rang_down = max(f_int):-int_param:min(f_int);

    % loop for bifurcation processing 2 (sweeping down)
    for i = 1:length(f_rang_down)

        % f value in step
        f = f_rang_down(i);

        % struc variable physical parameters
        X_var.ksi    = ksi;
        X_var.chi    = chi;
        X_var.lambda = lambda;
        X_var.kappa  = kappa;
        X_var.beta   = beta;
        X_var.f      = f;
        X_var.Omega  = Omega;
        X_var.delta  = delta;
        X_var.phi    = phi;

        % period of a forcing cycle
        T = 2*pi/Omega;

        % number of forcing cycles
        Nf = 5000;

        % initial dimensionless time
        t0 = 0.0;

        % final dimensionless time
        t1 = t0 + Nf*T;

        % number of samples per forcing cycle
        Nsamp = 10;

        % time series sampling points
        tspan = t0:(T/Nsamp):t1;

        % cut transient regime
        Nss = round(length(tspan)*0.7);
        Npp = length(tspan);

        % compute solve differential equation
        [~,Y1] = piezomagbeam_asymmetric(X_var,x0,tspan);

        % vector of bifurcation point
        Bifurc_down(:,i,k) = Y1(Nss:Nsamp:Npp,3);

        % rewrite initial condition
        x0 = [Y1(end,1) Y1(end,2) Y1(end,3)];
    end
    end

    disp(' ')
    disp(' ... ')
    disp(' ---- bifurcation processing "finished" --- ')
    disp(' ... ')

    % Save datas 
    disp(' ')
    disp(' ... ')
    disp(' ---- saving --- ')
    disp(' ... ')

    bifurc_inf.Bifurc_up   = Bifurc_up;
    bifurc_inf.Bifurc_down = Bifurc_down;
    bifurc_inf.Omega_rang  = Omega_rang;
    bifurc_inf.f_rang_up   = f_rang_up;
    bifurc_inf.f_rang_down = f_rang_down;

    save(name_file,'bifurc_inf')

    disp(' ')
    disp(' ... ')
    disp(' ---- finished --- ')
    disp(' ... ')
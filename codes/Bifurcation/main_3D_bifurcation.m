% -----------------------------------------------------------------
%  main_3D_bifurcation.m
% -----------------------------------------------------------------
%  This is the main file for a program which trace bifurcation 
%  diagram for the piezo-magneto-elastic beam%  
% -----------------------------------------------------------------
%  programmer: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Oct 20, 2020
% -----------------------------------------------------------------
%% Processing
clc
clear
close all


% program header
% -----------------------------------------------------------

disp(' ---------------------------------------------------')
disp('            Bifurcation Diagram of PMEH             ')
disp('                    (plot 3D)                       ')
disp('                                                    ')
disp(' by                                                 ')
disp(' João Pedro Norenberg (Unesp)                       ')
disp(' jpcvalese@gmail.com                                ')
disp(' ---------------------------------------------------')

% -----------------------------------------------------------

ksi    = 0.01;    % mechanical damping ratio
chi    = 0.05;    % dimensionless piezoeletric coupling term (mechanical)
lambda = 0.05;    % dimensionless time constant reciprocal
kappa  = 0.5;     % dimensionless piezoeletric coupling term (eletrical)
beta   = 2.5;     % nonlinear electromechanical coupling term

% struct variable physical parameters
var_input.X_params.ksi    = ksi;
var_input.X_params.chi    = chi;
var_input.X_params.lambda = lambda; 
var_input.X_params.kappa  = kappa;
var_input.X_params.beta   = beta;

% min and max of frequency
var_input.Par1_rang.Omega_rang = [0.64 0.96];

% size frequency vector
var_input.N1_rang.N_omega = 5;

% min and max of amplitude
var_input.Par2_rang.f_int = [0.02 0.30];

% increment of amplitude
var_input.N2_rang.int_param = 0.001;

% file name
name_file = ['bifurc_3D_beta',...
                    num2str(beta*10)];

% compute bifurcation
[bifurc_inf] = bifurcation_3d(var_input,name_file);

%% Plot
% plotting bifurcation (graph 3-D)
name_file = 'teste';
fig = plot_3D_bifurcation(bifurc_inf,name_file);
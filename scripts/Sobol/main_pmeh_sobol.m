% -----------------------------------------------------------------
%  main_pmeh_sobol.m
%
%  This script is the main file for a program that performs
%  a global sensitivity analysis in the nonlinear dynamics
%  of a Piezo-Magneto-Elastic-Beam Harvesting via Sobol indices.
% ----------------------------------------------------------------- 
%         João Pedro C V Norenberg
%            jpcvalese@gmail.com
%
%  last update: Dez 18, 2019
% -----------------------------------------------------------------
%  Reference:
%  
%  main_galloping_sobol.m
%  A global sensitivity analysis in the nonlinear dynamics
%  of a Galloping Oscillator via Sobol indices.
%
%  By: Americo Barbosa da Cunha Junior
clc
clearvars
close all

% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Piezo-Magneto-Elastic-Beam Harvesting              ')
disp(' (global sensitivity analysis)                      ')
disp('                                                    ')
disp(' by                                                 ')
disp(' João Pedro C V Norenberg                           ')
disp(' jpcvalese@gmail.com                                ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------


% start RNG and fix the statistical seed
% -----------------------------------------------------------
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream);
% -----------------------------------------------------------


% initialize UQLab
% -----------------------------------------------------------
uqlab
% -----------------------------------------------------------


% simulation information
% -----------------------------------------------------------
disp(' '); 
disp(' piezo.magneto.harvesting_sobol ');
disp(' ');

% -----------------------------------------------------------


% physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining model nominal parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');
   ksi1    = 0.01;       % damping ratio
   chi1    = 0.05;       % piezoelectric coupling (mechanical)
   lambda1 = 0.05;       % reciprocal time constant
   kappa1  = 0.5;        % piezoelectric coupling (electrical)
   f1      = 0.147;      % amplitude of external force
   Omega1  = 0.8;        % frequency of external force
% mathematical model
% -----------------------------------------------------------
disp(' '); 
disp(' --- defining mathematical model --- ');
disp(' ');
disp('    ... ');
disp(' ');

ModelOpts.mFile        = 'piezo_magneto_harvesting';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);
% -----------------------------------------------------------


% probabilistic input model
% -----------------------------------------------------------
disp(' ');
disp(' --- probabilistic input model --- ');
disp(' ');
disp('    ... ');
disp(' ');

% dispersion with respect to the nominal value
delta = 0.2;

InputOpts.Marginals(1).Name = 'ksi';
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [1-delta 1+delta]*ksi1;

InputOpts.Marginals(2).Name = 'chi';
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [1-delta 1+delta]*chi1;

InputOpts.Marginals(3).Name = 'lambda';
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [1-delta 1+delta]*lambda1;

InputOpts.Marginals(4).Name = 'kappa';
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [1-delta 1+delta]*kappa1;

InputOpts.Marginals(5).Name = 'f';
InputOpts.Marginals(5).Type = 'Uniform';
InputOpts.Marginals(5).Parameters = [1-delta 1+delta]*f1;

InputOpts.Marginals(6).Name = 'Omega';
InputOpts.Marginals(6).Type = 'Uniform';
InputOpts.Marginals(6).Parameters = [1-delta 1+delta]*Omega1;

myInput = uq_createInput(InputOpts);
%% -----------------------------------------------------------

% MC-based Sobol indices
% -----------------------------------------------------------

Order        = 1;
SampleSize   = 10;
repBootstrap = 5;
alpha        = 0.05;

SobolAnalysisMC = sobol_mc(Order,SampleSize,repBootstrap,alpha);

file_name = ['Sobol_MC_Nsamp',num2str(SampleSize),'_ord',num2str(Order),'_f',num2str(f1*1e3),...
                '_O',num2str(Omega1*1e1),'.mat'];
save(file_name,'SobolAnalysisMC')

plot_sobol(SobolAnalysisMC,1,'MC')
% -----------------------------------------------------------
%% -----------------------------------------------------------
% PCE-based Sobol indices
% -----------------------------------------------------------
degreePCE = 3;
Nsamples  = 100;
order     = 2;
SobolAnalysisPCE = sobol_pce(myInput,myModel,degreePCE,Nsamples,order);

file_name = ['Sobol_PCE_Nsamp',num2str(Nsamples),'_ord',num2str(order),'_f',num2str(f1*1e3),...
                '_O',num2str(Omega1*1e1),'.mat'];
save(file_name,'SobolAnalysisPCE')

plot_sobol(SobolAnalysisPCE,1,'PCE')
% -----------------------------------------------------------
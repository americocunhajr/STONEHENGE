% -----------------------------------------------------------------
%  sobol_mc.m
%
%  This function computes the Sobol Indecies to the piezo-magneto-
%  elastic beam based on MC method.
% ----------------------------------------------------------------- 
%  programmer: João Pedro C V Norenberg
%              jpcvalese@gmail.com
%
%  last update: Oct 21, 2019
% -----------------------------------------------------------------

function mySobolAnalysisMC = sobol_mc(Order,SampleSize,repBootstrap,alpha)

% MC-based Sobol indices
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- MC-based Sobol indices --- ');
disp(' ');
disp('    ... ');
disp(' ');

SobolOpts.Type             = 'Sensitivity';
SobolOpts.Method           = 'Sobol';
SobolOpts.Sobol.Order      = Order; 
SobolOpts.Sobol.SampleSize = SampleSize;

% Bootstrap Method
SobolOpts.Bootstrap.Replications = repBootstrap;    
SobolOpts.Bootstrap.Alpha        = alpha;

mySobolAnalysisMC = uq_createAnalysis(SobolOpts);

toc
% -----------------------------------------------------------
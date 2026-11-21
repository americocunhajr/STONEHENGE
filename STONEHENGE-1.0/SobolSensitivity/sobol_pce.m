% -----------------------------------------------------------------
%  sobol_pce.m
%
%  This function computes the Sobol Indecies to the piezo-magneto-
%  elastic beam based on PCE method.
% ----------------------------------------------------------------- 
%  programmer: João Pedro C V Norenberg
%              jpcvalese@gmail.com
%
%  last update: Oct 21, 2019
% -----------------------------------------------------------------

function mySobolAnalysisPCE = sobol_pce(myInput,myModel,degreePCE,Nsamples,order)
% PCE-based Sobol indices
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- PCE-based Sobol indices --- ');
disp(' ');
disp('    ... ');
disp(' ');

PCEOpts.Type               = 'Metamodel';
PCEOpts.MetaType           = 'PCE';
PCEOpts.Method             = 'LARS';
PCEOpts.Input              = myInput;
PCEOpts.FullModel          = myModel;
PCEOpts.Degree             = degreePCE;
PCEOpts.ExpDesign.NSamples = Nsamples;

myPCE = uq_createModel(PCEOpts);

PCESobol.Type        = 'Sensitivity';
PCESobol.Method      = 'Sobol';
PCESobol.Sobol.Order = order;

mySobolAnalysisPCE = uq_createAnalysis(PCESobol);
toc
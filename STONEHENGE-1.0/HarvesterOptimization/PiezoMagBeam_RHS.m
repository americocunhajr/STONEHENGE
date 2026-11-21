% -----------------------------------------------------------------
%  piezomagbeam.m
% -----------------------------------------------------------------
%  This function defines the right hand side of the follwoing
%  nonlinear system of ordinary differential equations
%
%    dy1/dt = y2
%    dy2/dt = -2*ksi.y2 + 0.5*y1(1-y1^2) + chi.y3 + f.cos(Omega.t)
%    dy3/dt = - lambda.y3 - kappa.y2,
%
%  that models the dynamics of a piezo-magneto-elastic beam.
%  
%  Reference:
%  A. Erturk, J. Hoffmann, and D. J. Inman
%  A piezomagnetoelastic structure for broadband vibration
%  energy harvesting
%  Applied Physics Letters
%  vol. 94 pp. 254102, 2009
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunha@uerj.br
%
%  last update: October 23, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function ydot = PiezoMagBeam_RHS(t,y,param)

% preallocate memory for Ydot vector
ydot = zeros(3,1);

% physical parameters
f      = param.f;
Omega  = param.Omega;
ksi    = param.ksi;
chi    = param.chi;
lambda = param.lambda;
kappa  = param.kappa;

% state space system of equations
ydot(1) = y(2);
ydot(2) = -2*ksi*y(2) + 0.5*y(1)*(1.0-y(1)^2) + chi*y(3) + f*cos(Omega*t);
ydot(3) = -lambda*y(3) - kappa*y(2);

end
% -----------------------------------------------------------------

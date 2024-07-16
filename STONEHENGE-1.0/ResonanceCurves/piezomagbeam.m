

% -----------------------------------------------------------------
%  piezomagbeam.m
%
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
%  
%  A. Erturk, J. Hoffmann, and D. J. Inman
%  A piezomagnetoelastic structure for broadband vibration
%  energy harvesting
%  Applied Physics Letters
%  vol. 94 pp. 254102, 2009
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Dec 8, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function ydot = piezomagbeam(t,y,phys_param)

% preallocate memory for Ydot vector
ydot = zeros(3,1);

% physical parameters
ksi    = phys_param(1);
chi    = phys_param(2);
f      = phys_param(3);
Omega  = phys_param(4);
lambda = phys_param(5);
kappa  = phys_param(6);
%x0     = phys_param(7);
%xdot0  = phys_param(8);
%v0     = phys_param(9);

% state space system of equations
ydot(1) = y(2);
ydot(2) = -2*ksi*y(2) + 0.5*y(1)*(1.0-y(1)^2) + chi*y(3) + f*cos(Omega*t);
ydot(3) = -lambda*y(3) - kappa*y(2);

end
% -----------------------------------------------------------------

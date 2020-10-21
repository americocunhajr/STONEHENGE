% -----------------------------------------------------------------
%  piezomagbeam_power.m
%
%  This function computes the power of a piezo-magneto-elastic
%  beam, which is defined as
%
%                  -- time(end)
%                 |
%  power := 1/T   | lambda Qvolt(t)^2 dt
%                 |
%               -- t = time(1)
%
%  where T = time(end) - time(1).
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
function [power,mean_power] = piezomagbeam_power(time,Qvolt,phys_param)

    % physical parameters
    %ksi    = phys_param(1);
    %chi    = phys_param(2);
    %f      = phys_param(3);
    %Omega  = phys_param(4);
    lambda = phys_param(5);
    %kappa  = phys_param(6);
    %x0     = phys_param(7);
    %xdot0  = phys_param(8);
    %v0     = phys_param(9);

    % temporal interval of analysis
    T = time(end)-time(1);

    % output power
    power = lambda*Qvolt.^2;
    mean_power = (1/T)*trapz(time,power);

end
% -----------------------------------------------------------------

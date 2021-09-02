% -----------------------------------------------------------------
%  BistableHarvesterPower.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% -----------------------------------------------------------------
%  This function computes the mean output power recovered from
%  a bistable piezo-magneto-elastic energy harvester, defined by
%
%                       -- time(end)
%                      |
%  power_mean := 1/T   | lambda volt(t)^2 dt
%                      |
%                    -- t = time(1)
%
%  where T = time(end) - time(1), volt is the output voltage,
%  and lambda is the reciprocal of circuit time constant.
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [power,power_mean] = HarvesterPower(time,volt,lambda)

% time interval of analysis
T = time(end)-time(1);

% output power
power = lambda*volt.^2;

% mean output power
power_mean = (1/T)*trapz(time,power);

end
% -----------------------------------------------------------------
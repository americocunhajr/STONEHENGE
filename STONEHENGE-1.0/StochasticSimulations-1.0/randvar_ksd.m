
% -----------------------------------------------------------------
%  randvar_ksd.m
%
%  This functions computes the kernel smooth density estimation
%  of a random process distribution given its numerical series.
%
%  input:
%  data   - (Ns x Ndt) data matrix
%  numpts - number of points
%
%  output:
%  data_ksd  - (numpts x Ndt) frequency counts matrix
%  data_supp - (numpts x Ndt) bins locations matrix
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 6, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [data_ksd,data_supp] = randvar_ksd(data,numpts)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    % compute data matrix dimensions
	[Ns,Ndt] = size(data);
    
    % check arguments
%    if isinteger(numpts) ~= 1
%    	error('numpts must be integer.')
%    end
    
    % preallocate memory for bins matrix
	data_supp = zeros(numpts,Ndt);
    
    % preallocate memory for frequency matrix
    data_ksd = zeros(numpts,Ndt);
    
    % loop over time instants
    for n=1:Ndt
        
        % maximum of data(:,n)
        data_max = max(data(:,n));
    
        % minimum of data(:,n)
        data_min = min(data(:,n));
        
        % points to be evaluated
        data_pts = linspace(data_min,data_max,numpts);
    
        % compute KSD estimation for data(:,n) distribution
        [data_ksd(:,n),data_supp(:,n)] = ksdensity(data(:,n),data_pts);
    end

return
% -----------------------------------------------------------------

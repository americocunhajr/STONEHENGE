
% -----------------------------------------------------------------
%  randvar_pdf.m
%
%  This functions computes the probability density funtion
%  of a random process given its numerical series.
%
%  input:
%  data     - (Ns x Ndt) data matrix 
%  numbins  - number of bins
%
%  output:
%  bins  - (numbins x Ndt) bins locations matrix (row vector)
%  freq  - (numbins x Ndt) frequency counts matrix (row vector)
%  area  -       (1 x Ndt) area under the histogram (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jun 23, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [bins,freq,area] = randvar_pdf(data,numbins)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    % compute data matrix dimensions
	[Ns,Ndt] = size(data);
    
    % check arguments
%    if isinteger(numbins) ~= 1
%    	error('numbins must be integer.')
%    end
    
	if numbins >= Ns
    	error('numbins must be less than Ns.')
	end
    
    % preallocate memory for bins matrix
	bins = zeros(numbins,Ndt);
    
    % preallocate memory for frequency matrix
    freq = zeros(numbins,Ndt);
    
    % preallocate memory for area vector
    area = zeros(1,Ndt);
    
    % loop over time instants
    for n=1:Ndt
    
        % maximum of data(:,n)
        data_max = max(data(:,n));
    
        % minimum of data(:,n)
        data_min = min(data(:,n));
        
        % check for potential errors
        if data_max == data_min
            
            if mod(numbins,2) == 0
                
                data_min = data_min - numbins/2;
                data_max = data_max + numbins/2 - 1;
                
            else
                
                data_min = data_min - (numbins-1)/2;
                data_max = data_max + (numbins-1)/2;
            end
        end
        
        % compute bin width
        binwidth = (data_max-data_min)/(numbins-1);
    
        % compute bins locations vector
        bins(:,n) = (data_min:binwidth:data_max);

        % compute frequency counts vector
        freq(:,n) = histc(data(:,n),bins(:,n));
        
        % normalize frequency counts matrix
        freq(:,n) = freq(:,n)/(Ns*binwidth);
        
        % compute area under the histogram
        area(1,n) = binwidth*sum(freq(:,n));
    end

return
% -----------------------------------------------------------------

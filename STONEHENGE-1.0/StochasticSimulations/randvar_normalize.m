
% -----------------------------------------------------------------
%  randvar_normalize.m
%
%  This function normalizes the numerical series of a
%  given random process.
%
%  input:
%  X     - (Ns x Ndt) random process numerical series
%
%  output:
%  Xnorm - normalized random variable
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Jun 23, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function Xnorm = randvar_normalize(X)

    % check number of arguments
    if nargin < 1
        error('Too few inputs.')
    elseif nargin > 1
        error('Too many inputs.')
    end
    
    % compute data matrix dimensions
	[Ns,Ndt] = size(X);
    
    % preallocate memory for normalized matrix
    Xnorm = zeros(Ns,Ndt);
    
    % loop over time instants
    for n=1:Ndt
        
        % normalize the numerical series
        Xnorm(:,n) = (X(:,n) - mean(X(:,n))) / std(X(:,n));
    end

return
% -----------------------------------------------------------------

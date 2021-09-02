% -----------------------------------------------------------------
%  fredholm_expcorr_eig.m
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: May 10, 2016
% ----------------------------------------------------------------- 
%  This is function computes the eigenpairs of Fredholm integral
%  operator, with kernel K(x1,x2) = exp(-|x1-x2|/a), where a > 0
%  is a correlation length. This kernel is the covariance function
%  of an underlying real-valued random field {S(x), x \in D}, with
%  domain D = [-b,b]. The corresponding Fredholm equation is
%  
%     --
%    |
%    | K(x1,x2) phi_n(x2) dx2 = lambda_n phi_n(x1), for x1 in D
%    |
%  -- D
%  
%  
%  References:
%  R. Ghanem and P. Spanos
%  Stochastic Finite Element Method: A Spectral Approch
%  Dover Publications, 2003
%  Pages 28-33
%  
%  D. Xiu
%  Numerical Methods for Stochastic Computations
%  Princeton University Press, 2009
%  Pages 47-49
%  
%  
%  Input:
%  b        - domain upper limit
%  a        - autocovariance correlation length
%  Neig     - number of eigenpairs to be computed
%  Nx       - number of mesh points for domain discretization
%  tol      - numerical tolerance for root find  (optional)
%  max_iter - maximum of iteration for root find (optional)
%
%  Output:
%  lambda  - ( 1 x Neig) vector with eigenvalues
%  phi     - (Nx x Neig) matrix with eigenfunction (in columns)
%  xmesh   - ( 1 x Nx  ) domain discretization mesh
%  dphid x - (Nx x Neig) matrix with eigenfunction derivatives (in columns)
%  omega   - ( 1 x Neig) vector with charac. eq. solutions
%  iroot   - charac. eq. solutions counter
%  ising   - charac. eq. singularities counter
%  iter    - iterations counter
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [lambda,phi,xmesh,dphidx,omega,iroot,ising,iter] = ...
                 fredholm_expcorr_eig(b,a,Neig,Nx,tol,max_iter)

if nargin < 4
    error('Too few inputs.')
elseif nargin > 6
    error('Too many inputs.')
elseif nargin == 4
    % tolerance
    tol = 1.0e-7;
    % maximum of iterations
    max_iter = 1e6;
elseif nargin == 5
    % maximum of iterations
    max_iter = 1e6;
end

% charact. eq. root lower bound
v1_odd  = 0.0;
w1_even = 0.0;

% charact. eq.root upper bound
v2_odd  = 0.0;
w2_even = 0.0;

% discretization refinement parameter
M = 100;

% mesh size
% (using correlation length inverse as characteristic value)
dw = (1/a)/M;

% root counter
iroot = 0;

% singularity counter
ising = 0;

% iteration counter
iter = 0;


% preallocate memory for characteristic equation solutions
omega = zeros(1,Neig);

% preallocate memory for eigenvalues
lambda = zeros(1,Neig);

% preallocate memory for eigenfunctions
phi = zeros(Nx,Neig);

% preallocate memory for eigenfunctions derivatives
dphidx = zeros(Nx,Neig);

% domain discretization mesh
xmesh = linspace(-b,b,Nx);


% autocorrelation function eigenpairs
while iter < max_iter && iroot < Neig
    
    % update iteration counter
    iter = iter + 1;
    
    % check if next root index is odd or even
    if mod(iroot+1,2) ~= 0
        
        % update upper and lower (odd) bounds
        v1_odd = v2_odd;
        v2_odd = v1_odd + dw;
        
        % evaluate charact. eq. at the (odd) bounds
        f1 = 1-a*v1_odd*tan(b*v1_odd);
        f2 = 1-a*v2_odd*tan(b*v2_odd);
    else
        
        % update upper and lower (even) bounds 
        w1_even = w2_even;
        w2_even = w1_even + dw;
        
        % evaluate charact. eq. at the (even) bounds
        f1 = a*w1_even+tan(b*w1_even);
        f2 = a*w2_even+tan(b*w2_even);
    end
        
	% check the necessary condition to be a root
    if f1*f2 < 0
    	
    	% check if next root index is odd or even
        if mod(iroot+1,2) ~= 0
            
            % look for a root (odd case)
            [v_odd,fv_odd] = fzero(@(v_odd) 1-a*v_odd*tan(b*v_odd),[v1_odd,v2_odd]);
    
            % check if v_odd is a root or a singularity
            if abs(fv_odd) < tol
            
				% update root index
                iroot = iroot + 1;
                
                % save charact. eq. solutions
                omega(1,iroot) = v_odd;
                
                % corresponding eigenvalue
                lambda(1,iroot) = (2*a)./(1+a^2*v_odd^2);
                
                % normalization coefficient
                NC = 1./sqrt(b+sin(2*v_odd*b)/(2*v_odd));
                
                % corresponding eigenfunction
                    phi(:,iroot) =        cos(v_odd*xmesh)*NC;
                 dphidx(:,iroot) = -v_odd*sin(v_odd*xmesh)*NC;
            else
                
				% update singulatity counter
                ising = ising + 1;
            end
        
        else
            
            % look for a root (even case)
            [w_even,fw_even] = fzero(@(w_even) a*w_even+tan(b*w_even),[w1_even,w2_even]);
    
            % check if w_even is a root or a singularity
            if abs(fw_even) < tol
                
				% update root index
                iroot = iroot + 1;
				
                % save charact. eq. solutions
                omega(1,iroot) = w_even;
                    
                % corresponding eigenvalue
                lambda(1,iroot) = (2*a)./(1+a^2*w_even^2);
                
                % normalization coefficient
                NC = 1./sqrt(b-sin(2*w_even*b)/(2*w_even));
               
                % corresponding eigenfunction
                    phi(:,iroot) =        sin(w_even*xmesh)*NC;
                 dphidx(:,iroot) = w_even*cos(w_even*xmesh)*NC;
            else
                
				% update singulatity counter
                ising = ising + 1;
            end
        end
    end
end

return
% -----------------------------------------------------------

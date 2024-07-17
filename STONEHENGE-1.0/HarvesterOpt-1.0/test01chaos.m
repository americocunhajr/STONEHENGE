% -----------------------------------------------------------------
%  test01chaos.m
% -----------------------------------------------------------------
%  This functions applies Gottwald-Melbourne 0-1 test for chaos
%  to a given time series. Result is near to 0 for non-chaotic 
%  data and near 1 for chaotic data.
%  
%  Adapted from Z1TEST.m by Paul Matthews, July 2009
%  http://arxiv.org/pdf/0906.1418v1
%  https://goo.gl/vWSqXs
%  
%  References:
%  Georg A. Gottwald and Ian Melbourne
%  The 0-1 Test for Chaos: A review
%  In Chaos Detection and Predictability,
%  Springer Lecture Notes in Physics 915 
%  Editors: C. Skokos, G.A. Gottwald and J. Laskar, 2016
%  http://www.springer.com/gp/book/9783662484081
%  
%  D. Bernardini and G. Litak
%  An overview of 0-1 test for chaos
%  J Braz. Soc. Mech. Sci. Eng. (2016) 38: 1433.
%  doi:10.1007/s40430-015-0453-y
%
%  input:
%  x        - (N x 1) observable time series
%  cmin     - lower bound for test parameter c
%  cmax     - upper bound for test parameter c
%  Nc       - number of test repeats
%  plotflag - diagnostic plots flag (1 - turn on / 0 - turn off)
%  OSflag   - oversampling msg flag (1 - turn on / 0 - turn off)
%
%  output:
%  kmedian - ( 1 x  1) classifier (correlation vector mediam)
%  kcorr   - (Nc x  1) correlation vector
%  c       - (Nc x  1) parameter c vector
%  p       - ( N x Nc) generalized dispalcement p
%  q       - ( N x Nc) generalized dispalcement q
%  OS_ID1  - oversampling indicator 1
%  OS_ID2  - oversampling indicator 2
% -----------------------------------------------------------------
%  programmer: Americo Cunha
%              americo.cunha@uerj.br
%
%  last update: October 23, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [kmedian,kcorr,c,p,q,OS_ID1,OS_ID2] = ...
                test01chaos(x,cmin,cmax,Nc,OSflag,plotflag)

% check number of arguments
if nargin < 1
    error('Too few inputs.')
elseif nargin > 6
    error('Too many inputs.')
elseif nargin == 1
        cmin = 0.0;
        cmax = 2*pi;
          Nc = 100;
      OSflag = 1;
    plotflag = 0;
elseif nargin == 2
    error('[cmin,cmax] interval should be prescribed.')
elseif nargin == 3
          Nc = 100;
      OSflag = 1;
    plotflag = 0;
elseif nargin == 4
      OSflag = 1;
    plotflag = 0;
elseif nargin == 5
    plotflag = 0;
end

% check if c interval is valid
if cmax <= cmin || cmin < 0.0 || cmax > 2*pi
    error('[cmin,cmax] must be a subset of [0,2*pi] interval.')
end

% compute the time series dimenions
s = size(x);

% convert time series to a column vector (if necessary)
if s(1) == 1
    x = x';
end

% number of points in time series
N = length(x);

% number of points in short time series
N0 = round(N/10);

% check time series size
if N0 < 1
    %error('time series is too short')
    N0 = N;
end

% index vetor for original time series
% (N x 1) vector
j = (1:N)';

% random values for c in [cmin,cmax] \subset [0,2*pi]
% (Nc x 1) vector
c = cmin + rand(Nc,1)*(cmax-cmin);

% square of observable time average
AVG2 = mean(x)^2;

% reduced time series to represent the line y = x
% (N0 x 1) vector
bisec_line = (1:N0)';

% preallocate memory for generalized coordinates
p = zeros(N,Nc);
q = zeros(N,Nc);

% preallocate memory for mean square displacement
% (N0 x 1) vector
MSD = zeros(N0,1);

% prealloctate memory for correlation vector (classifier)
% (Nc x 1) vector
kcorr = zeros(Nc,1);

% loop over samples of parameter c
for it=1:Nc
    
   % transform time series from (x,xdot) to (p,q) space
   p(:,it) = cumsum(x.*cos(j*c(it,1)));
   q(:,it) = cumsum(x.*sin(j*c(it,1)));
   
   % mean square displacement minus oscilatory part
   for n=1:N0;
      
      
      MSD(n,1) = mean((p(n+1:N)-p(1:N-n)).^2  + ...
                      (q(n+1:N)-q(1:N-n)).^2) - ...
                 AVG2*(1-cos(n*c(it,1)))/(1-cos(c(it,1)));
   end
   
   % correlation coefficient
   kcorr(it,1) = corr(bisec_line,MSD);
end

% correlation vector mediam (classifier)
kmedian = abs(median(kcorr));

% two crude attempts to check for oversampling
OS_ID1 = (max(x)-min(x))/mean(abs(diff(x)));
OS_ID2 = median(kcorr(c<mean(c))) - median(kcorr(c>mean(c)));

% oversampling message
if OSflag == 1
    
    if OS_ID1 > 10 || OS_ID2 > 0.5
        disp('Warning: data is probably oversampled.')
        disp('Use coarser sampling or reduce the maximum value of c.')
    end
end

% useful diagnostic plots
if plotflag == 1
    
    figure(1);
    plot(c,abs(kcorr),'*');
    xlabel('c');ylabel('Kc');
    
    figure(2)
    plot(bisec_line,MSD);
    xlabel('N');ylabel('MSD');
    
    figure(3)
    plot(p(:,1),q(:,1),'o-r','MarkerEdgeColor','b');
    xlabel('p');ylabel('q');
end

return
% -----------------------------------------------------------------

% Z1TEST  Gottwald-Melbourne 0-1 test for chaos
% 
% Z1TEST(X) is the result of the 0-1 test applied to the vector X.
% Result is near to 0 for non-chaotic data and near 1 for chaotic data.
% see  http://arxiv.org/pdf/0906.1418v1.  Paul Matthews, July 2009.

function [kmedian, N_over]=z1test(x)

s = size(x); if s(2)==1; x=x'; end
N = length(x); j=1:N;
t = 1:round(N/10);
M = zeros(1,round(N/10));
ns = 100; % Number of samples
c = rand(1,ns)*pi;  % ns random c values in [0,pi]
kcorr = zeros(1,ns);

N_over = 0;

for its=1:ns
   p=cumsum(x.*cos(j*c(its)));q=cumsum(x.*sin(j*c(its)));
   for n=1:round(N/10)
      M(n)=mean( (p(n+1:N)-p(1:N-n)).^2 + (q(n+1:N)-q(1:N-n)).^2 )- ...
           mean(x)^2*(1-cos(n*c(its)))/(1-cos(c(its)));
   end
   kcorr(its)=corr(t',M');
end

%  plot(c,kcorr,'*');xlabel('c');ylabel('k'); % useful diagnostic plots
%  plot(t,M);xlabel('t');ylabel('M')

% Two crude attempts to check for oversampling:
 if (max(x)-min(x) )/mean(abs(diff(x))) > 10 || ...
        median(kcorr(c<mean(c))) - median(kcorr(c>mean(c))) > 0.5
     %disp('Warning: data is probably oversampled.')
     %disp('Use coarser sampling or reduce the maximum value of c.')
     N_over = 1;
 end

kmedian=abs(median(kcorr));
end

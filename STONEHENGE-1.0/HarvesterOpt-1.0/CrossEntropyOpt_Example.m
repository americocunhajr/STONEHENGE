
% -----------------------------------------------------------------
%  One-dimensional optimization with the Cross-entropy method
%  
%  Example from:
%  Reuven Y. Rubinstein, Dirk P. Kroese
%  The Cross-Entropy Method: A Unified Approach to
%  Combinatorial Optimization, Monte-Carlo Simulation,
%  and Machine Learning
%  Springer-Verlag, 2004.
% -----------------------------------------------------------------

clc
clear
close all

S = inline('exp(-(x-2).^2) + 0.8*exp(-(x+2).^2)');

mu = -6;
sigma = 10;
Nel = 10;
N = 100;
eps = 1E-8;
tic
t=0;

while sigma > eps
t = t+1;
x = mu + sigma*randn(N,1);
SX = S(x);
sortSX = sortrows([x SX],2);
Xel = sortSX((N - Nel + 1):N,1);
mu = mean(Xel);
sigma = std(Xel);
fprintf('%g %6.9f %6.9f %6.9f \n', t, S(mu),mu, sigma)
end

mu
toc

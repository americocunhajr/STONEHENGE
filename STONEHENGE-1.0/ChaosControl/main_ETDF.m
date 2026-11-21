% -----------------------------------------------------------------
%  main_ETDF.m
% ----------------------------------------------------------------- 
%  This is the main file for a program for the piezo-magneto-elastic 
%  beam with control chaos by the Extended time-delayed feedback method.
% ----------------------------------------------------------------- 
%  programmers: 
%        José Geraldo Telles Ribeiro (telles@eng.uerj.br)
%        
% -----------------------------------------------------------------


clear all
close all


% Simulation parameters
tf=10000;
dt=0.001;
t=0:dt:tf-dt;

% Input signal
A=1.3*9.81;
u=(0.03*A)*1.4;
w=0.8;

% initial conditions
x0=[-1.63 0.78 0 0];
x0=[1 0 0];

% Duffing Parameters
si=0.01;
k=0.51;
chi=0.51;
lam=0.04;

% Digital controller parameters
T=0.01;
K=-1.3;
R=0.95;
% Run Simulink
sim Harvest_DuffingCL


% Plotting
figure('color',[1 1 1]);
subplot(2,2,1)
plot(x(length(t)/2:length(t),1),x(length(t)/2:length(t),2),'linewidth',2.0);grid;
ax = gca;
ax.FontSize = 16;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];
ax.LineWidth = [2];
xlabel('Displacement');ylabel('Velocity');
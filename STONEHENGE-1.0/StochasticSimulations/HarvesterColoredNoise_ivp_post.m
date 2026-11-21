% -----------------------------------------------------------------
%  HarvesterColoredNoise_ivp_post.m
% -----------------------------------------------------------------
%  programmer: Americo Cunha Jr
%              americo.cunhajr@gmail.com
%
%  last update: Sep 7, 2020
% -----------------------------------------------------------------
%  Script for post processing of the simulation data from
%  HarvesterColoredNoise_ivp_main.m
% -----------------------------------------------------------------


% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% plot external force
% ...........................................................
gtitle = ' ';
leg1 = 'force';
leg2 = 'noise';
xlab   = ' time';
ylab   = ' external force';
xmin   = t0;
xmax   = t1;
%ymin   = -0.15;
%ymax   = 0.15;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__force_vs_time'];
gname2  = [num2str(case_name),'__force2_vs_time'];
flag   = 'eps';
fig01  = graph_type2(time_noise,force_harmonic,...
                     time_noise,noise,...
                     gtitle,leg1,leg2,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname2,flag);

fig0b  = graph_type1(time_noise,force_total,gtitle,....
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig0a);
% ...........................................................



% plot displacement
% ...........................................................
xlab   = 'time';
ylab   = 'displacement';
gtitle = ' ';
xmin   = t0;
xmax   = t1;
%ymin   = -2.0;
%ymax   =  2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__disp_vs_time'];
flag   = 'eps';
fig1a  = graph_type1(time,Qdisp,gtitle,xlab,ylab,...
                     xmin,xmax,ymin,ymax,gname,flag);
%close(fig1a);
% ...........................................................



% plot velocity
% ...........................................................
xlab   = 'time';
ylab   = 'velocity';
gtitle = '';
xmin   = t0;
xmax   = t1;
%ymin   = -2.0;
%ymax   =  2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__velo_vs_time'];
flag   = 'eps';
fig1b  = graph_type1(time,Qvelo,gtitle,xlab,ylab,...
                     xmin,xmax,ymin,ymax,gname,flag);
%close(fig1b);
% ...........................................................


% plot voltage
% ...........................................................
xlab   = 'time';
ylab   = 'voltage';
gtitle = '';
xmin   = t0;
xmax   = t1;
%ymin   = -2.0;
%ymax   =  2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__volt_vs_time'];
flag   = 'eps';
fig1c  = graph_type1(time,Qvolt,gtitle,xlab,ylab,...
                     xmin,xmax,ymin,ymax,gname,flag);
%close(fig1c);
% ...........................................................



% plot displacement (steady state)
% ...........................................................
% xlab   = 'time';
% ylab   = 'displacement';
% gtitle = 'steady state';
% xmin   = time(Nss);
% xmax   = time(Ndt);
% %ymin   = -2.0;
% %ymax   =  2.0;
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__disp_vs_time_ss'];
% flag   = 'eps';
% fig2a  = graph_type1(time(Nss:Ndt),Qdisp(Nss:Ndt),gtitle,xlab,ylab,...
%                      xmin,xmax,ymin,ymax,gname,flag);
%close(fig1a);
% ...........................................................



% plot velocity (steady state)
% ...........................................................
% xlab   = 'time';
% ylab   = 'velocity';
% gtitle = 'steady state';
% xmin   = time(Nss);
% xmax   = time(Ndt);
% %ymin   = -2.0;
% %ymax   =  2.0;
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__velo_vs_time_ss'];
% flag   = 'eps';
% fig2b  = graph_type1(time(Nss:Ndt),Qvelo(Nss:Ndt),gtitle,xlab,ylab,...
%                      xmin,xmax,ymin,ymax,gname,flag);
%close(fig2b);
% ...........................................................


% plot voltage (steady state)
% ...........................................................
% xlab   = 'time';
% ylab   = 'voltage';
% gtitle = 'steady state';
% xmin   = time(Nss);
% xmax   = time(Ndt);
% %ymin   = -2.0;
% %ymax   =  2.0;
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__volt_vs_time_ss'];
% flag   = 'eps';
% fig2c  = graph_type1(time(Nss:Ndt),Qvolt(Nss:Ndt),gtitle,xlab,ylab,...
%                      xmin,xmax,ymin,ymax,gname,flag);
%close(fig2c);
% ...........................................................



% plot power time series
% ...........................................................
xlab   = 'time';
ylab   = 'power';
leg1   = 'power';
leg2   = 'mean';
gtitle = '';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__power_vs_time'];
flag   = 'eps';
fig2d  = graph_type2(time,power,...
                     time,power_mean*ones(1,Ndt),...
                     gtitle,leg1,leg2,xlab,ylab,...
                     xmin,xmax,ymin,ymax,gname,flag);
%close(fig2d);
% ...........................................................





% plot phase trajectoty disp vs velo
% ...........................................................
% xlab   = 'displacement';
% ylab   = 'velocity';
% gtitle = 'phase space trajectory';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
% %xmin   = -2.0;
% %xmax   =  2.0;
% %ymin   = -2.0;
% %ymax   =  2.0;
% gname  = [num2str(case_name),'__phasespace_disp_vs_velo'];
% flag   = 'eps';
% fig3a  = plot_2d_trajectory(Qdisp(Nss:Ndt),Qvelo(Nss:Ndt),...
%                             gtitle,xlab,ylab,...
%                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig3a);
% ...........................................................


% plot phase trajectoty disp vs volt
% ...........................................................
% xlab   = 'displacement';
% ylab   = 'voltage';
% gtitle = 'phase space trajectory';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
% %xmin   = -2.0;
% %xmax   =  2.0;
% %ymin   = -1.0;
% %ymax   =  1.0;
% gname  = [num2str(case_name),'__phasespace_disp_vs_volt'];
% flag   = 'eps';
% fig3b  = plot_2d_trajectory(Qdisp(Nss:Ndt),Qvolt(Nss:Ndt),...
%                             gtitle,xlab,ylab,...
%                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig3b);
% ...........................................................


% plot phase trajectoty velo vs volt
% ...........................................................
% xlab   = 'velocity';
% ylab   = 'voltage';
% gtitle = 'phase space trajectory';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
% %xmin   = -2.0;
% %xmax   =  2.0;
% %ymin   = -1.0;
% %ymax   =  1.0;
% gname  = [num2str(case_name),'__phasespace_velo_vs_volt'];
% flag   = 'eps';
% fig3c  = plot_2d_trajectory(Qvelo(Nss:Ndt),Qvolt(Nss:Ndt),...
%                             gtitle,xlab,ylab,...
%                             xmin,xmax,ymin,ymax,gname,flag);
%close(fig3c);
% ...........................................................


% plot displacement histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' displacement';
ylab   = ' probability density';
%xmin   = -5.0;
%xmax   =  5.0;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
%ymax   = 2.5;
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_disp'];
flag   = 'eps';
fig4a  = graph_bar_curve1(Qdisp_bins,Qdisp_freq,...
                          Qdisp_supp,Qdisp_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4a);
% -----------------------------------------------------------


% plot velocity histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' velocity';
ylab   = ' probability density';
%xmin   = -5.0;
%xmax   =  5.0;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
%ymax   = 2.5;
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_velo'];
flag   = 'eps';
fig4b  = graph_bar_curve1(Qvelo_bins,Qvelo_freq,...
                          Qvelo_supp,Qvelo_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4b);
% -----------------------------------------------------------



% plot voltage histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' voltage';
ylab   = ' probability density';
%xmin   = -5.0;
%xmax   =  5.0;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 0.0;
%ymax   = 2.5;
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_volt'];
flag   = 'eps';
fig4c  = graph_bar_curve1(Qvolt_bins,Qvolt_freq,...
                          Qvolt_supp,Qvolt_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4c);
% -----------------------------------------------------------

% plot power histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' power';
ylab   = ' probability density';
xmin   = 0.0;
%xmax   =  5.0;
xmax   = 'auto';
ymin   = 0.0;
%ymax   = 2.5;
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_power'];
flag   = 'eps';
fig4c  = graph_bar_curve1(power_bins,power_freq,...
                          power_supp,power_ksd,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4c);
% -----------------------------------------------------------


% plot dynamical system atractor
% ...........................................................
Njump = round(0.1*(1/dt));

xlab   = 'displacement';
ylab   = 'velocity';
zlab   = 'voltage';
gtitle = 'dynamical system atractor';
xmin   = min(Qdisp(1:Ndt));
xmax   = max(Qdisp(1:Ndt));
ymin   = min(Qvelo(1:Ndt));
ymax   = max(Qvelo(1:Ndt));
zmin   = min(Qvolt(1:Ndt));
zmax   = max(Qvolt(1:Ndt));
gname  = [num2str(case_name),'__atractor'];
flag   = 'eps';
fig5   = plot_3d_atractor(Qdisp(1:Njump:Ndt),...
                          Qvelo(1:Njump:Ndt),...
                          Qvolt(1:Njump:Ndt),...
                           time(1:Njump:Ndt),...
                          gtitle,xlab,ylab,zlab,...
                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig5);
% ...........................................................


toc
% -----------------------------------------------------------
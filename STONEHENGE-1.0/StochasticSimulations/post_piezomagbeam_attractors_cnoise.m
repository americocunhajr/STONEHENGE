
% -----------------------------------------------------------------
%  Script for post processing of the simulation data from
%  main_piezomagbeam_attractors_cnoise.m
% -----------------------------------------------------------------


% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% plot displacement
% ...........................................................
xlab   = 'time';
ylab   = 'displacement';
gtitle = ' ';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
xmin   = t0;
xmax   = t1;
%ymin   = -2.0;
%ymax   =  2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__disp_vs_time'];
flag   = 'eps';
fig1a  = graph_type3(time1,Qdisp1,...
                     time2,Qdisp2,...
                     time3,Qdisp3,...
                     gtitle,leg1,leg2,leg3,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1a);
% ...........................................................



% plot velocity
% ...........................................................
xlab   = 'time';
ylab   = 'velocity';
gtitle = '';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
xmin   = t0;
xmax   = t1;
%ymin   = -2.0;
%ymax   =  2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__velo_vs_time'];
flag   = 'eps';
fig1b  = graph_type3(time1,Qvelo1,...
                     time2,Qvelo2,...
                     time3,Qvelo3,...
                     gtitle,leg1,leg2,leg3,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1b);
% ...........................................................


% plot voltage
% ...........................................................
xlab   = 'time';
ylab   = 'voltage';
gtitle = '';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
xmin   = t0;
xmax   = t1;
%ymin   = -2.0;
%ymax   =  2.0;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__volt_vs_time'];
flag   = 'eps';
fig1c  = graph_type3(time1,Qvolt1,...
                     time2,Qvolt2,...
                     time3,Qvolt3,...
                     gtitle,leg1,leg2,leg3,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1c);
% ...........................................................


% plot power time series
% ...........................................................
xlab   = 'time';
ylab   = 'power';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
gtitle = '';
xmin   = t0;
xmax   = t1;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__power_vs_time'];
flag   = 'eps';
fig2d  = graph_type3(time1,power1,...
                     time2,power2,...
                     time3,power3,...
                     gtitle,leg1,leg2,leg3,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig2d);
% ...........................................................



% plot phase trajectoty disp vs velo
% ...........................................................
xlab   = 'displacement';
ylab   = 'velocity';
gtitle = 'phase space trajectory';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
%xmin   = -2.0;
%xmax   =  2.0;
%ymin   = -2.0;
%ymax   =  2.0;
gname  = [num2str(case_name),'__phasespace_disp_vs_velo'];
flag   = 'eps';
%fig3a  = plot_2d_trajectory(Qdisp(Nss1:Ndt),Qvelo(Nss1:Ndt),...
%                            gtitle,xlab,ylab,...
%                            xmin,xmax,ymin,ymax,gname,flag);
%close(fig3a);
% ...........................................................


% plot phase trajectoty disp vs volt
% ...........................................................
xlab   = 'displacement';
ylab   = 'voltage';
gtitle = 'phase space trajectory';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
%xmin   = -2.0;
%xmax   =  2.0;
%ymin   = -1.0;
%ymax   =  1.0;
gname  = [num2str(case_name),'__phasespace_disp_vs_volt'];
flag   = 'eps';
%fig3b  = plot_2d_trajectory(Qdisp(Nss1:Ndt),Qvolt(Nss1:Ndt),...
%                           gtitle,xlab,ylab,...
%                            xmin,xmax,ymin,ymax,gname,flag);
%close(fig3b);
% ...........................................................


% plot phase trajectoty velo vs volt
% ...........................................................
xlab   = 'velocity';
ylab   = 'voltage';
gtitle = 'phase space trajectory';
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
%xmin   = -2.0;
%xmax   =  2.0;
%ymin   = -1.0;
%ymax   =  1.0;
gname  = [num2str(case_name),'__phasespace_velo_vs_volt'];
flag   = 'eps';
%fig3c  = plot_2d_trajectory(Qvelo(Nss1:Ndt),Qvolt(Nss1:Ndt),...
%                            gtitle,xlab,ylab,...
%                            xmin,xmax,ymin,ymax,gname,flag);
%close(fig3c);
% ...........................................................


% plot displacement histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' displacement';
ylab   = ' probability density';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
%xmin   = -5.0;
%xmax   =  5.0;
%ymin   = 0.0;
%ymax   = 2.5;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_disp'];
flag   = 'eps';
fig4a  = graph_bar_curve3(Qdisp1_bins,Qdisp1_freq,...
                          Qdisp2_bins,Qdisp2_freq,...
                          Qdisp3_bins,Qdisp3_freq,...
                          Qdisp1_supp,Qdisp1_ksd,...
                          Qdisp2_supp,Qdisp2_ksd,...
                          Qdisp3_supp,Qdisp3_ksd,...
                          leg1,leg2,leg3,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4a);
% -----------------------------------------------------------


% plot velocity histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' velocity';
ylab   = ' probability density';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
%xmin   = -5.0;
%xmax   =  5.0;
%ymin   = 0.0;
%ymax   = 2.5;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_velo'];
flag   = 'eps';
fig4b  = graph_bar_curve3(Qvelo1_bins,Qvelo1_freq,...
                          Qvelo2_bins,Qvelo2_freq,...
                          Qvelo3_bins,Qvelo3_freq,...
                          Qvelo1_supp,Qvelo1_ksd,...
                          Qvelo2_supp,Qvelo2_ksd,...
                          Qvelo3_supp,Qvelo3_ksd,...
                          leg1,leg2,leg3,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4b);
% -----------------------------------------------------------



% plot voltage histogram
% -----------------------------------------------------------
gtitle = ' ';
xlab   = ' voltage';
ylab   = ' probability density';
leg1   = 'non-zero disp';
leg2   = 'non-zero velo';
leg3   = 'non-zero volt';
%xmin   = -5.0;
%xmax   =  5.0;
%ymin   = 0.0;
%ymax   = 2.5;
xmin   = 'auto';
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__pdf_volt'];
flag   = 'eps';
fig4c  = graph_bar_curve3(Qvolt1_bins,Qvolt1_freq,...
                          Qvolt2_bins,Qvolt2_freq,...
                          Qvolt3_bins,Qvolt3_freq,...
                          Qvolt1_supp,Qvolt1_ksd,...
                          Qvolt2_supp,Qvolt2_ksd,...
                          Qvolt3_supp,Qvolt3_ksd,...
                          leg1,leg2,leg3,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4c);
% -----------------------------------------------------------


% plot dynamical system atractor
% ...........................................................
Njump = 1;

xlab   = 'displacement';
ylab   = 'velocity';
zlab   = 'voltage';
gtitle = 'dynamical system atractor';
xmin   = min(Qdisp1(Nss1:Ndt));
xmax   = max(Qdisp1(Nss1:Ndt));
ymin   = min(Qvelo1(Nss1:Ndt));
ymax   = max(Qvelo1(Nss1:Ndt));
zmin   = min(Qvolt1(Nss1:Ndt));
zmax   = max(Qvolt1(Nss1:Ndt));
%xmin   = -2.0;
%xmax   =  2.0;
%ymin   = -2.0;
%ymax   =  2.0;
%zmin   = -1.0;
%zmax   =  1.0;
gname  = [num2str(case_name),'__atractor'];
flag   = 'eps';
%fig5   = plot_3d_atractor(Qdisp1(Nss1:Njump:Ndt),...
%                          Qvelo1(Nss1:Njump:Ndt),...
%                          Qvolt1(Nss1:Njump:Ndt),...
%                           time(Nss2:Njump:Ndt),...
%                          gtitle,xlab,ylab,zlab,...
%                          xmin,xmax,ymin,ymax,zmin,zmax,gname,flag);
%close(fig5);
% ...........................................................


toc
% -----------------------------------------------------------
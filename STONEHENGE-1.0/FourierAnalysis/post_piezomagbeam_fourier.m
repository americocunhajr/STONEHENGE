
% -----------------------------------------------------------------
%  Script for post processing of the simulation data from
%  main_piezomagbeam_fourier.m
% -----------------------------------------------------------------



% -----------------------------------------------------------
tic

disp(' ');
disp(' --- post processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% plot displacement PSD
% ...........................................................
gtitle = ' displacement';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 0.0;
xmax   = freq_max;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_disp'];
flag   = 'eps';
fig4a  = graph_type2_psd(freq,10*log10(psd_Qdisp),...
                         freq,10*log10(psd_Qdisp_sg),...
                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4a);
% ...........................................................


% plot velocity PSD
% ...........................................................
gtitle = ' velocity';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 0.0;
xmax   = freq_max;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_velo'];
flag   = 'eps';
fig4b  = graph_type2_psd(freq,10*log10(psd_Qvelo),...
                         freq,10*log10(psd_Qvelo_sg),...
                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4b);
% ...........................................................


% plot voltage PSD
% ...........................................................
gtitle = ' voltage';
xlab   = ' frequency (Hz)';
ylab   = ' power spectral density (dB/Hz)';
xmin   = 0.0;
xmax   = freq_max;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__psd_volt'];
flag   = 'eps';
fig4c  = graph_type2_psd(freq,10*log10(psd_Qvolt),...
                         freq,10*log10(psd_Qvolt_sg),...
                         gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4c);
% ...........................................................


toc
% -----------------------------------------------------------

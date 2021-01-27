% -----------------------------------------------------------------
%  ploting.m
% -----------------------------------------------------------------
%  This is the main file for a program which plot the Basins of 
%  attraction attraction for the piezo-magneto-elastic beam
% -----------------------------------------------------------------
%  programmer: 
%        João Peterson (ligier.peterson@gmail.com)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Jan, 2021
% -----------------------------------------------------------------

%clc
%clear
%close all

xdim = 1200;
ydim = 1200;

buffer = zeros(xdim,ydim);

a = 2;
f = 19+(a-1)*32;

IC_name = "IC_[-3,3]_[-3,3]_" + xdim + "X" + ydim + "_f" + f + "_O80";
AT_name = "dat_files/AT_[-3,3]_[-3,3]_" + xdim + "X" + ydim + "_f" + f + "_O80";

filename_IC = IC_name + ".dat";
IC = importdata(filename_IC);

for i=1:xdim
   buffer(i,:) = IC(end-i+1,:); 
end

x0_min = -3;
x0_max = 3;
N_x0 = xdim;

xdot0_min = -3;
xdot0_max = 3;
N_xdot0 = ydim;

x0 = linspace(x0_min, x0_max, N_x0);
xdot0 = linspace(xdot0_min, xdot0_max, N_xdot0);

x      = x0;
y      = xdot0;
F      = buffer;
gtitle = '';
xlab   = 'initial displacement';
ylab   = 'initial velocity';
xmin   = x0_min;
xmax   = x0_max;
ymin   = xdot0_min;
ymax   = xdot0_max;
gname  = char(IC_name);
flag   = 'png';

fig = graph_contourf_pnt(x,y,F,gtitle,xlab,ylab,...
                                 xmin,xmax,ymin,ymax,gname,flag);
close(fig);

filename_AT = AT_name + ".dat";
AT = importdata(filename_AT);
color = [ 0.0 1.0 0.0;  % green
          1.0 0.0 0.0;  % red
          0.0 0.0 1.0;  % blue
          0.0 1.0 1.0;  % cyan
          1.0 0.0 1.0;  % magenta
          1.0 1.0 0.0;  % yellow
          %0.5 1.0 1.0;  %
          0.0 0.0 1.0;  % blue
          %1.0 0.5 1.0;  %
          1.0 0.0 0.0;  % red
          %1.0 1.0 0.5   %
          0.0 1.0 0.0;  % green
          %0.3 0.5 0.8   %
          0.0 1.0 0.0;  % green
          ];

for b=1
   % if b==10
        plot(AT(0.8*end:end,3*b-2), AT(0.8*end:end,3*b-1), 'Color', color(b,:));
   % else
    %    plot(AT(:,3*b-2), AT(:,3*b-1), 'Color', color(b,:));
    %end
    xlabel('displacement', 'FontSize', 16, 'FontName', 'Helvetica');
    ylabel('velocity', 'FontSize', 16, 'FontName', 'Helvetica');
    %        xlim([-3,3]);
    %        ylim([-3,3]);
    set(gcf,'color','white');
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);

    hold on 
end
hold off

gname = char(AT_name);
saveas(gcf, gname, 'png')
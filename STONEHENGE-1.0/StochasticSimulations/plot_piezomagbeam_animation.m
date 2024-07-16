
% -----------------------------------------------------------------
%  plot_piezomagbeam_animation.m
%
%  This functions plots an animation for the nonlinear dynamics
%  of an piezo-magneto-elastic beam.
%
%  input:
%  time        - (1 x Ndt) time vector
%  y1          - (1 x Ndt) beam tip displacement
%  B1          - trailer base left arm length (m)
%  B2          - trailer base right arm length (m)
%  L1          - trailer vertical arm length (m)
%  L1          - tower arm length (m)
%  Dtire       - tires diamenter (m)
%  vtrans_km_h - velocity of translation (km/h)
%  vtitle      - video title
%  legend      - legend text
%  xmin        - x axis minimum value
%  xmax        - x axis maximum value
%  ymin        - y axis minimum value
%  ymax        - y axis maximum value
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Dec 11, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = plot_piezomagbeam_animation(time,Qdisp,vtitle,legend,...
                                           xmin,xmax,ymin,ymax)
    
    % check number of arguments
    if nargin < 8
        error('Too few inputs.')
    elseif nargin > 8
        error('Too many inputs.')
    end

    % check arguments
    if length(time) ~= length(Qdisp)
        error('vectors time and Qdisp must have same length')
    end
    
    % convert to row vector (if necessary)
    if find( size(time) == max(size(time)) ) < 2
        time=time';
    end
    
    if find( size(Qdisp) == max(size(Qdisp)) ) < 2
        Qdisp=Qdisp';
    end
	
    
    % number of time steps
    Ndt = length(time);
    
    % domain dimensions
    xlength = xmax-xmin;
    ylength = ymax-ymin;
    
    % rigid base dimensions
    rb_length = 0.50*xlength;
    rb_heigth = 0.80*ylength;
    rb_thick  = 0.10*rb_length;
    
    % beam dimension
    beam_heigth = 0.75*rb_heigth;
    
    % mesh for beam mode shape
    zmesh = beam_heigth*linspace(0.0,1.0,10);
    
    % beam mode shape
    b1 = (0.5*pi + 0.3042)*300e-3;
    phi1 = sinh(b1*zmesh) - sin(b1*zmesh) - ...
           ((sinh(b1)+sin(b1))/(cosh(b1)+cos(b1)))*(cosh(b1*zmesh)-cos(b1*zmesh));
    
    % coordinates of the rigid base vertex
    xv1 = xmin + 0.25*xlength;
    yv1 = ymin;
    xv2 = xv1;
    yv2 = yv1 + rb_heigth;
    xv3 = xv2 + rb_length;
    yv3 = yv2;
    xv4 = xv3;
    yv4 = yv3 - rb_thick;
    xv5 = xv4 - rb_length + rb_thick;
    yv5 = yv4;
    xv6 = xv5;
    yv6 = yv5 - rb_heigth + 2*rb_thick;
    xv7 = xv6 + rb_length - rb_thick;
    yv7 = yv6;
    xv8 = xv7;
    yv8 = yv7 - rb_thick;
    
    xrb = [xv1 xv2 xv3 xv4 xv5 xv6 xv7 xv8];
    yrb = [yv1 yv2 yv3 yv4 yv5 yv6 yv7 yv8];

    % dashed vertical line coordinates
    xvline1 = xv1 + 0.55*rb_length;
    xvline2 = xv1 + 0.55*rb_length;
    yvline1 = yv4;
    yvline2 = yv4 - beam_heigth;
    
    % magnets coordinates
    xm1 = xv1 + 0.25*rb_length;
    ym1 = yv6;

    xm2 = xv1 + 0.75*rb_length;
    ym2 = yv6;
    
    % resistor coordinates
    xr1 = xvline1 - 0.1*rb_length;
    yr1 = 1.2*yv2;
    xr2 = xvline1 - 0.15*rb_length;
    yr2 = 0.7*yv4;
    
    % coordinates of the time counter
    xtime = 0.8*xmin;
    ytime = 0.85*ymax;
    
    
    % loop to construct the video
    for n=1:Ndt
        
        % initialize video frame
        fig = figure(1000);
        
        % define frame properties
        set(gcf,'color','white');
        set(gca,'Box','on');
        set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
        set(gca,'XMinorTick','off','YMinorTick','off');
        set(gca,'XGrid','off','YGrid','off');
        set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
        set(gca,'FontName','Helvetica');
        set(gca,'FontSize',16);
        
        sp1 = subplot(2,1,1);
        
        % draw rigid base bottom part
        fh01 = fill(xrb,yrb,[0.95,0.95,0.95]);
        
        hold on
        
        % draw dashed vertical line
        fh02 = plot([xvline1 xvline2],[yvline1 yvline2],'k--');
        
        % draw magnets
        fh03 = rectangle('Position',[xm1,ym1,0.10*rb_length,0.10*rb_length],...
                         'FaceColor',[0,0,0]);
        fh04 = rectangle('Position',[xm2,ym2,0.10*rb_length,0.10*rb_length],...
                         'FaceColor',[0,0,0]);
        
        % draw resistor
        fh05 = rectangle('Position',[xr1,yr1,0.2*rb_length,0.05*rb_length],...
                         'FaceColor',[0,0,0]);
        fh06 = rectangle('Position',[xr2,yr2,0.3*rb_length,0.40*rb_length]);
        
        % draw beam
        fh07 = plot(xvline1+phi1*Qdisp(n),yvline1 - zmesh,'b-','LineWidth',6);
        
        % display time counter
        text(xtime,ytime,['time = ',num2str(time(n),'%.1f')],...
                          'Color','k',...
                          'FontName','Helvetica',...
                          'FontSize',14);

                      
        % define video title
        title(vtitle,'FontSize',16,'FontName','Helvetica');
        
        % no axis trick
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        
        % define video frame limits
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        
        hold off
        
        sp2 = subplot(2,1,2);
        
        % plot displacement time series
        fh08 = plot(time,Qdisp,'b-',time(n),Qdisp(n),'*r');
        set(gca,'FontSize',16);
        
        % define video title
        xlabel('time','FontSize',16,'FontName','Helvetica')
        ylabel('displacement','FontSize',16,'FontName','Helvetica')
        
        
        
        set(sp1, 'Position', [0.15, 0.35, 0.7, 0.55])
        set(sp2, 'Position', [0.15, 0.10, 0.7, 0.2])
        
        
        % pause exibition
        if n == 1
            %pause
        end

    end

end
% -----------------------------------------------------------------

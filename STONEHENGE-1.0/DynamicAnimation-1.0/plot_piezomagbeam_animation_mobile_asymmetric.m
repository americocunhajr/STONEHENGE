% -----------------------------------------------------------------
%  plot_piezomagbeam_animation_mobile_asymmetric.m
%
%  This functions plots an animation for the nonlinear dynamics
%  of an asymmetric piezo-magneto-elastic beam in non-inertial frame
%  of reference.
%
%  input:
%  time        - (1 x Ndt) time vector
%  vtitle      - video title
%  magflag     - magnet on/off flag
%  xmin        - x axis minimum value
%  xmax        - x axis maximum value
%  ymin        - y axis minimum value
%  ymax        - y axis maximum value
% ----------------------------------------------------------------- 
%  programmers: 
%         Americo Cunha Jr (americo.cunhajr@gmail.com)
%         Joao Pedro Norenberg (jpcvalese@gmail.com)  
%
%  last update: March 30, 2021
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function plot_piezomagbeam_animation_mobile_asymmetric(series, f  ,Omega,...
                                                 dim,beta,  IC,phi,delta)
warning('off')
    %check number of arguments
    if nargin < 8
        error('Too few inputs.')
    elseif nargin > 8
        error('Too many inputs.')
    end
    
    % struct input
    time  = series.time;
    Qdisp = series.Disp;
    Qvolt = series.Volt;
    xmin  = dim.xmin;
    xmax  = dim.xmax;
    ymin  = dim.ymin;
    ymax  = dim.ymax;

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
    xlength = xmax - xmin;
    ylength = ymax - ymin;
    
    % rigid base dimensions
    rb_length = 0.50 * xlength;
    rb_heigth = 0.80 * ylength;
    rb_thick  = 0.10 * rb_length;
    
    % matrix rotation
    T = [cos(phi.*pi/180) -sin(phi.*pi/180)
         sin(phi.*pi/180)  cos(phi.*pi/180)];
    
    % beam dimension
    beam_heigth = 0.75 * rb_heigth;
    
    % mesh for beam mode shape
    zmesh = beam_heigth*linspace(0.0,1.0,10);
    
    % beam mode shape
    b1 = (0.5*pi + 0.3042)*300e-3;
    phi1 = sinh(b1*zmesh) - sin(b1*zmesh) - ...
           ((sinh(b1)+sin(b1))/(cosh(b1)+cos(b1)))*(cosh(b1*zmesh)-cos(b1*zmesh));
    
    % initial coordinates of the rigid base vertex
    xv1_0 = xmin + 0.25*xlength;
    
    % dashed vertical line coordinates
    xvline1_0 = xv1_0 + 0.55*rb_length;
    
    % coordinates of the time counter
    xtime = 0.8*xmin;
    ytime = 0.9*ymax;
    
    % file save name
    v_name = ['f',num2str(f*1e3),'_O',num2str(Omega*1e1)  ,  '_ic', ...
              num2str(IC(1)),num2str(IC(2)),num2str(IC(3)),'_beta', ...
              num2str(beta*10),'_delta',num2str(delta*100),'_phi' , ...
              num2str(phi),'_mobile'];
    
    % open video file
    myVideo           = VideoWriter(v_name,'MPEG-4');
    myVideo.Quality   = 100;  
    myVideo.FrameRate = 40;
    open(myVideo) 
    
    Njump = 50;
    
    % loop to construct the video
     for n=1:Njump:Ndt 
        
        % initialize video frame
        fig = figure(1000);
        set(fig,'Position',[500 50 580 660]);
        
        % define frame properties
        set(gcf,'color','white');
        set(gca,'Box','on');
        set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
        set(gca,'XMinorTick','off','YMinorTick','off');
        set(gca,'XGrid','off','YGrid','off');
        set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
        set(gca,'FontName','Helvetica');
        set(gca,'FontSize',16);
        
        
        % cossine of Omega t
		cos_Omega_t = 1; %f*cos(Omega*time(n));
		
		% coordinates of the rigid base vertex
		xv1 = xv1_0*cos_Omega_t;
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
		
        desl_v = 0;%0.5;
        
	    xrb = [xv1 xv2 xv3 xv4 xv5 xv6 xv7 xv8]-desl_v*ones(1,8);
    	yrb = [yv1 yv2 yv3 yv4 yv5 yv6 yv7 yv8];
        A1 = [xrb; yrb];
        
        % rotation base
        base_rot = T*A1;
        
        
    	% dashed vertical line coordinates
        xvline1 = xv1 + 0.55*rb_length - desl_v;
    	xvline2 = xvline1;
        yvline1 = yv4;
        yvline2 = yv4 - beam_heigth;
        
        
        sp1 = subplot(3,1,1);
        
        % draw rigid base bottom part
        fh01 = fill(base_rot(1,:),base_rot(2,:),[0.8,0.8,0.8]);
       
        hold on
        
        fh13 = rectangle('Position' ,[xv1_0-0.2,-1.11,1.4,0.01], ...
                         'EdgeColor',[0.25 0.25 0.25]         , ...
                         'FaceColor',[0.25 0.25 0.25] );
       
        % draw dashed vertical line
        P1 = T*[xvline1;yvline1];
        P2 = T*[xvline2;yvline2];
        fh02 = plot([P1(1,1) P2(1,1)],[P1(2,1) P2(2,1)],'k--');
        
        % magnets coordinates
        xm1 = xv1 + 0.25*rb_length+0.1;
        ym1 = yv6;

        % draw magnets
        fh03 = patch('XData',[xm1,xm1,xm1+0.10*rb_length,xm1+0.10*rb_length],...
                     'YData',[ym1,ym1+0.10*rb_length,ym1+0.10*rb_length,ym1], ...
                     'EdgeColor',  [0.9100 0.4100 0.1700]               , ...
                     'FaceColor',  [0.9100 0.4100 0.1700] );
                 
        rotate(fh03,[0 0 1],phi)

        xm2 = xv1 + 0.75*rb_length;
        ym2 = yv6;
        
        fh04 = patch('XData',[xm2,xm2,xm2+0.10*rb_length,xm2+0.10*rb_length],...
                     'YData',[ym2,ym2+0.10*rb_length,ym2+0.10*rb_length,ym2], ...
                     'EdgeColor',  [0.9100 0.4100 0.1700]               , ...
                     'FaceColor',  [0.9100 0.4100 0.1700] );
        rotate(fh04,[0 0 1],phi)
		
	    % resistor coordinates
    	xr1 = xvline1 - 0.1*rb_length+0.08;
	    yr1 = 1.2*yv2;
    	xr2 = xvline1 - 0.15*rb_length+0.08;
	    yr2 = 0.7*yv4;
        
        % draw resistor
        %fh06 = rectangle('Position',[xr2,yr2,0.3*rb_length,0.40*rb_length]);
        fh05 = patch('XData',[xr2,xr2,xr2+0.3*rb_length,xr2+0.3*rb_length],...
                     'YData',[yr2,yr2+0.40*rb_length,yr2+0.40*rb_length,yr2],...
                     'FaceAlpha',0);
        rotate(fh05,[0 0 1],phi)
		
        fh06 = patch('XData',[xr1,xr1,xr1+0.2*rb_length,xr1+0.2*rb_length],...
                     'YData',[yr1,yr1+0.05*rb_length,yr1+0.05*rb_length,yr1],...
                     'FaceColor',[0,1,0]);               
        rotate(fh06,[0 0 1],phi)
        
        
        % draw pzt       
        coord_pzt = T*[xvline1+phi1(1:4)*Qdisp(n); yvline1 - zmesh(1:4)]; 
        fh10 = plot(coord_pzt(1,:), coord_pzt(2,:),...
                    '-','color',[0,0,1],'LineWidth',9.5);
        
        % draw beam
        coord_beam = T*[xvline1+phi1*Qdisp(n); yvline1 - zmesh];
        fh07 = plot(coord_beam(1,:),coord_beam(2,:), ...
                    'color' , [0.45,0.45,0.45]            , ...
                    'LineWidth', 6);

                
        % nonlinear coupling
        if beta ~= 0
            if n > 1
                delete(fh11);
            end
            Cnc = T*[0.57-cos_Omega_t*0.17 0.61-cos_Omega_t*0.17; 0.81 0.835];
            Cnc = [0.51 0.56; 0.82 0.850];
            fh11 = annotation('arrow',Cnc(1,:), Cnc(2,:), ...
                              'color' , [1,0,0], 'HeadLength' , 5   , ...
                              'HeadWidth'  , 7 ,'linewidth',2);
                          
            text(-0.3+xmin,1.1*ymax,['$\;$ asymmetric bistable energy harvester (nonlinear coupling)'], ...
                                        'Color','k','FontName','Helvetica','FontSize',14                 , ...
                                        'interpreter','latex');
        else
            text(-0.3+xmin,1.1*ymax                                                 , ...
                 '$\;\;$ asymmetric bistable energy harvester (linear coupling)' , ...
                  'Color', 'k', 'FontName','Helvetica', 'FontSize',14               , ...
                  'interpreter','latex');
        end
        
        % display time counter
        text(-0.23+xmin,ytime,[' time         = ' , num2str(time(n),'%.1f') ,...
                               ' $\;\;f$      = ' , num2str(f,'%.3f')       ,...
                               ' $\;\;\Omega$ = ' , num2str(Omega,'%.1f')   ,...
                               ' $\;\;\beta$  = ' , num2str(beta,'%.2f')    ,...
                               ' $\;\;IC$     = (', num2str(IC(1)) ,','     ,...
                                                    num2str(IC(2)) ,','     ,...
                                                    num2str(IC(3)) ,')' ]   ,...
                               'Color' , 'k', 'FontName','Helvetica'        ,...
                               'FontSize',14,'interpreter','latex');
           
        text(-0.705*xmin,-0.7*ymin,['$\phi$ = ', num2str(phi), '$^{\circ}$'],...
                                    'Color' , 'k', 'FontName','Helvetica'   ,...
                                    'FontSize',14,'interpreter','latex');
                                
        text(-0.715*xmin,-0.5*ymin,['$\delta$ = ', num2str(delta)]          ,...
                                    'Color' , 'k', 'FontName','Helvetica'   ,...
                                    'FontSize',14,'interpreter','latex');
                           
        text(-0.6*xmin,0.4*ymin,['non-inertial', newline , ...
                                 '  frame of ' , newline , ...
                                 ' reference']);
                             
        % no axis trick
        set(gca,'xtick'  , [])
        set(gca,'ytick'  , [])
        set(gca,'visible','off')
        
        % define video frame limits
        xlim([xmin xmax]);
        ylim([ymin-0.2 ymax]);
        
        hold off
        
        sp2 = subplot(3,1,2);
        
        % plot displacement time series
        fh08 = plot(time,Qdisp,'b-',time(n),Qdisp(n),'or', ...
                    'MarkerFaceColor','r','MarkerSize',5);
        set(gca,'FontSize',16);
        
        % set axis
        ylabel('displacement','FontSize',14,'FontName','Helvetica')
        set(gca,'xtick',[])
        ylim([-2.4 2.4])
        
        sp3 = subplot(3,1,3);
        
        % plot voltage time series
        fh09 = plot(time,Qvolt,'b-',time(n),Qvolt(n),'or', ...
                    'MarkerFaceColor','r','MarkerSize',5);
        set(gca,'FontSize',16);
        
        % set axis
        xlabel('time','FontSize',14,'FontName','Helvetica')
        ylabel('voltage','FontSize',14,'FontName','Helvetica')
        volt_max = max(Qvolt);
        ylim([-3 3])
        
        % Reference (article and authors)
        text(-0.2,-0.7,'Global sensitivity analysis of (a)symmetric energy harvesters (2021)', ...
            'Units','normalized','FontName','Times New Roman','color',[.4 0 .4]           , ...
            'FontWeight', 'Bold')
        
        text(-0.2,-0.88,'J. P. Norenberg, A. Cunha Jr, S. da Silva, P. S. Varoto'         , ...
            'Units', 'normalized', 'FontName', 'Times New Roman', 'color',[.4 0 .4]       , ...
            'FontWeight', 'Bold')
        
        % position plots
        set(sp1, 'Position', [0.2, 0.55, 0.7, 0.40])
        set(sp2, 'Position', [0.2, 0.35, 0.7, 0.15])
        set(sp3, 'Position', [0.2, 0.17, 0.7, 0.15])
        
        
        pause(0.0001)
        
        % pause exibition
        if n == 1
            pause
        end
        
        F(n) = getframe(fig);
        writeVideo(myVideo, F(n));
     end
    myVideo.close()
 end
% -----------------------------------------------------------------

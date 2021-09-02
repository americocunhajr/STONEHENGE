
% -----------------------------------------------------------------
%  plot_3d_atractor.m
%
%  This functions plots a 3D atractor for a dynamical system
%  which time series are given.
%
%  input:
%  x1     - x data vector 1
%  y1     - y data vector 1
%  z1     - z data vector 1
%  time   - time vector
%  gtitle - graph title
%  xlab   - x axis label
%  ylab   - y axis label
%  ylab   - z axis label
%  xmin   - x axis minimum value
%  xmax   - x axis maximum value
%  ymin   - y axis minimum value
%  ymax   - y axis maximum value
%  zmin   - z axis minimum value
%  zmax   - z axis maximum value
%  gname  - graph name
%  flag   - output file format (optional)
%
%  output:
%  gname.eps - output file in eps format (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 10, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = plot_3d_atractor(x1,y1,z1,time,gtitle,xlab,ylab,zlab,...
                                xmin,xmax,ymin,ymax,zmin,zmax,gname,flag)
    
    % check number of arguments
    if nargin < 15
        error('Too few inputs.')
    elseif nargin > 16
        error('Too many inputs.')
    elseif nargin == 15
        flag = 'none';
    end

    % check arguments
    if length(x1) ~= length(y1) || ...
            length(x1) ~= length(z1) || ...
                    length(y1) ~= length(z1)
        error('vectors x1, y1, and z1 must be same length')
    end
    
    if length(x1) ~= length(time)
        error('x1 and time vectors must be same length')
    end
    
	global xdata ydata zdata tdata fh th
    
    xdata = x1;
    ydata = y1;
    zdata = z1;
    tdata = time;

    fig = figure('Name',gname,'NumberTitle','off');
    
    set(gcf,'color','white');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'zColor',[.3 .3 .3]);
    %set(gca,'Box','on');
    
    if ( strcmp(xmin,'auto') || strcmp(xmax,'auto') )
        xlim('auto');
    else
        xlim([xmin xmax]);
    end
    
    if ( strcmp(ymin,'auto') || strcmp(ymax,'auto') )
        ylim('auto');
    else
        ylim([ymin ymax]);
    end
    
    if ( strcmp(zmin,'auto') || strcmp(zmax,'auto') )
        zlim('auto');
    else
        zlim([zmin zmax]);
    end
    
    view(45, 45);

    hold on
    grid on

    % time display
    th = text(0.05,0.05,['time = ',num2str(time(1),'%.1f')],...
                         'Color','k','FontName',...
                         'Helvetica','FontSize',12,...
                         'Units','normalized');
    
    % plot initial data
    fh = plot3(x1(1),y1(1),z1(1));
                         
    set(fh,'LineWidth',1.0);
    set(fh,'MarkerSize',2.0);
    set(fh,'MarkerFaceColor','w');
    set(fh,'MarkerEdgeColor','k');

    xlabel(xlab,'FontSize',22,'FontName','Helvetica');
    ylabel(ylab,'FontSize',22,'FontName','Helvetica');
    zlabel(zlab,'FontSize',22,'FontName','Helvetica');
            
    title(gtitle,'FontSize',22,'FontName','Helvetica');

    % time series length
    Ndata = numel(time);
    
    % build timer
    TimerObj = timer('Period',0.001,... %period
                     'ExecutionMode','fixedRate',... %{singleShot,fixedRate,fixedSpacing,fixedDelay}
                     'BusyMode','drop',... %{drop, error, queue}
                     'TasksToExecute',Ndata-1,...
                     'StartDelay',0,...
                     'TimerFcn',@tcb,...
                     'StartFcn',[],...
                     'StopFcn',[],...
                     'ErrorFcn',[]);
      
    % start the timer
    start(TimerObj);
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        gname = [gname, '.eps'];
        %graph_fixPSlinestyle(gname,gname);
    end
    
end
% -----------------------------------------------------------------


% TimerFcn function
% -----------------------------------------------------------------
function tcb(src,evt)
      
    global xdata ydata zdata tdata fh th

    %What task are we on?  
    %Use this instead of for-loop variable ii
    taskEx = get(src,'TasksExecuted');

    % update the data.
    set(fh,'XData',xdata(1:taskEx),...
           'YData',ydata(1:taskEx),...
           'ZData',zdata(1:taskEx));
          
	% update time display
	set(th,'String',['time = ',num2str(tdata(taskEx),'%.1f')]);

	%force event queue flush
	drawnow;
end
% -----------------------------------------------------------------

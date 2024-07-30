
% -----------------------------------------------------------------
%  graph_type2.m
%
%  This functions plots a graph with two curves.
%
%  input:
%  x1     - x data vector 1
%  y1     - y data vector 1
%  x2     - x data vector 2
%  y2     - y data vector 2
%  gtitle - graph title
%  leg1   - legend 1
%  leg2   - legend 2
%  xlab   - x axis label
%  ylab   - y axis label
%  xmin   - x axis minimum value
%  xmax   - x axis maximum value
%  ymin   - y axis minimum value
%  ymax   - y axis maximum value
%  gname  - graph name
%  flag   - output file format (optional)
%
%  output:
%  gname.eps - output file in eps format (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Mar 21, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_type2(x1,y1,x2,y2,...
                           gtitle,leg1,leg2,...
                           xlab,ylab,xmin,xmax,ymin,ymax,gname,flag)
    
    % check number of arguments
    if nargin < 14
        error('Too few inputs.')
    elseif nargin > 15
        error('Too many inputs.')
    elseif nargin == 14
        flag = 'none';
    end

    % check arguments
    if length(x1) ~= length(y1)
        error('x1 and y1 vectors must be same length')
    end
    
    if length(x2) ~= length(y2)
        error('x2 and y2 vectors must be same length')
    end
        
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = plot(x1,y1,'-b');
    hold all
    fh2 = plot(x2,y2,'-r');
    set(gcf,'color','white');
    set(gca,'position',[0.2 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    %set(gca,'XTick',xmin:xmax);
    %set(gca,'YTick',ymin:ymax);
    %axis([xmin xmax ymin ymax]);
    
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
    
    set(fh1,'LineWidth',1.0);
    set(fh1,'MarkerSize',4.0);
    set(fh1,'MarkerFaceColor','w');
    set(fh1,'MarkerEdgeColor','b');
    set(fh2,'LineWidth',2.0);
    set(fh2,'MarkerSize',4.0);
    set(fh2,'MarkerFaceColor','w');
    set(fh2,'MarkerEdgeColor','r');
    leg = legend(leg1,leg2);
    set(leg,'FontSize',12);
    %set(leg,'interpreter', 'latex');
    labX = xlabel(xlab,'FontSize',16,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',16,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
    hold off
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        %gname = [gname, '.eps'];
        %graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------
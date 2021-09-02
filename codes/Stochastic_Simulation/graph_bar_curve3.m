
% -----------------------------------------------------------------
%  graph_bar_curve3.m
%
%  This functions plots a bar histogram and a curve with 
%  one data serie.
%
%  input:
%  bins1  - bins locations vector 1
%  freq1  - frequency counts vector 1
%  bins2  - bins locations vector 2
%  freq2  - frequency counts vector 2
%  bins3  - bins locations vector 3
%  freq3  - frequency counts vector 3
%  x1     - x data vector 1
%  y1     - y data vector 1
%  x2     - x data vector 2
%  y2     - y data vector 2
%  x3     - x data vector 3
%  y3     - y data vector 3
%  gtitle - graph title
%  leg1   - legend 1
%  leg2   - legend 2
%  leg3   - legend 3
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
%  last update: Sep 3, 2017
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_bar_curve3(bins1,freq1,bins2,freq2,bins3,freq3,...
                                x1,y1,x2,y2,x3,y3,gtitle,leg1,leg2,leg3,...
                                xlab,ylab,xmin,xmax,ymin,ymax,gname,flag)
	
    % check number of arguments
    if nargin < 23
        error('Too few inputs.')
    elseif nargin > 24
        error('Too many inputs.')
    elseif nargin == 23
        flag = 'none';
    end

    % check arguments
    if length(bins1) ~= length(freq1)
        error('bins1 and freq1 vectors must be same length')
    end
    
    if length(bins2) ~= length(freq2)
        error('bins2 and freq2 vectors must be same length')
    end
    
    if length(bins3) ~= length(freq3)
        error('bins3 and freq3 vectors must be same length')
    end
    
    if length(x1) ~= length(y1)
        error('x1 and y1 vectors must be same length')
    end
    
    if length(x2) ~= length(y2)
        error('x2 and y2 vectors must be same length')
    end
    
    if length(x3) ~= length(y3)
        error('x3 and y3 vectors must be same length')
    end
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    fh1 = bar(bins1,freq1,0.8);
    hold all
    fh2 = bar(bins2,freq2,0.8);
    fh3 = bar(bins3,freq3,0.8);
    fh4 = plot(x1,y1,'-b');
    fh5 = plot(x2,y2,'-m');
    fh6 = plot(x3,y3,'-k');
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
    
    set(fh1,'FaceColor','none');
    set(fh1,'EdgeColor','b');
    set(fh1,'LineStyle','-');
    set(fh2,'FaceColor','none');
    set(fh2,'EdgeColor','m');
    set(fh2,'LineStyle','-');
    set(fh3,'FaceColor','none');
    set(fh3,'EdgeColor','k');
    set(fh3,'LineStyle','-');
    set(fh4,'LineWidth',2.0);
    set(fh4,'MarkerSize',2.0);
    set(fh4,'MarkerFaceColor','w');
    set(fh4,'MarkerEdgeColor','k');
    set(fh5,'LineWidth',2.0);
    set(fh5,'MarkerSize',2.0);
    set(fh5,'MarkerFaceColor','w');
    set(fh5,'MarkerEdgeColor','k');
    set(fh6,'LineWidth',2.0);
    set(fh6,'MarkerSize',2.0);
    set(fh6,'MarkerFaceColor','w');
    set(fh6,'MarkerEdgeColor','k');
    leg = legend(leg1,leg2,leg3,3);
    labX = xlabel(xlab,'FontSize',18,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',18,'FontName','Helvetica');
    %set(Xlab,'interpreter','latex');
    %set(Ylab,'interpreter','latex');
    uistack(fh3,'top');
    uistack(fh2,'top');
    uistack(fh1,'top');
    
    hold off
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        %gname = [gname, '.eps'];
        %graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------

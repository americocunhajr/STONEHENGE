
% -----------------------------------------------------------------
%  graph_contourf_pnt.m
%
%  This functions plots the contour map of a scalar function
%  F: R^2 -> R.
%
%  input:
%  x      - x mesh vector
%  y      - y mesh vector
%  F      - scalar field
%  gtitle - graph title
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
%  last update: Oct 29, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_contourf_pnt(x,y,F,Xmax,Ymax,gtitle,xlab,ylab,...
                                  xmin,xmax,ymin,ymax,gname,flag)
                                    
    % check number of arguments
    if nargin < 13
        error('Too few inputs.')
    elseif nargin > 14
        error('Too many inputs.')
    elseif nargin == 13
        flag = 'none';
    end
    
    fig = figure('Name',gname,'NumberTitle','off');
    
    % generate 2D mesh grid
    [Xq,Yq] = meshgrid(x,y);
    fh = contourf(Xq,Yq,F);
    %colormap jet;
    %colormap cool
    %colormap pink
    %colormap hot
    colormap bone
    colorbar
    %scaxis([min(min(F)),max(max(F))])
    hold on
    fh2 = plot(Xmax,Ymax,'xb');
    hold off
    set(gcf,'color','white');
    set(gca,'position',[0.15 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','off');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'YDir','normal')
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',18);
    set(fh2,'LineWidth',2.0);
    set(fh2,'MarkerSize',10.0);
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
    labX = xlabel(xlab,'FontSize',18,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',18,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
	title(gtitle,'FontSize',20,'FontName','Helvetica');
    
    

    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        %gname = [gname, '.eps'];
        %graph_fixPSlinestyle(gname,gname);
    end

return
% -----------------------------------------------------------------
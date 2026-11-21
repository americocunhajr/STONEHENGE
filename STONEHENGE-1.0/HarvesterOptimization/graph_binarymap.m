% -----------------------------------------------------------------
%  graph_binarymap.m
% -----------------------------------------------------------------
%  This functions plots the contour map of a binary scalar 
%  function F: R^2 -> {0,1}.
%
%  input:
%  x      - x mesh vector
%  y      - y mesh vector
%  label0 - 0 value label
%  label1 - 1 value label
%  F      - scalar field
%  gtitle - graph title
%  xlab   - x axis label
%  ylab   - y axis label
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
%  programmer: Americo Cunha
%              americo.cunha@uerj.br
%
%  last update: March 31, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = graph_binarymap(x,y,F,Xmax,Ymax,gtitle,xlab,ylab,label0,label1,...
                              xmin,xmax,ymin,ymax,zmin,zmax,gname,flag)
                                    
    % check number of arguments
    if nargin < 17
        error('Too few inputs.')
    elseif nargin > 19
        error('Too many inputs.')
    elseif nargin == 17
        flag = 'none';
    end
    

    fig = figure('Name',gname,'NumberTitle','off');
    
    fh = imagesc(x,y,F);
    %colormap jet;
    %colormap bone
    colormap([1 1 1;0 0 0]) % black/white
    colorbar
    caxis([zmin zmax])
    colorbar('YTick',[0 1],'YTicklabel',{label0,label1})
    hold on
    fh2 = plot(Xmax,Ymax,'xr');
    hold off
    set(gcf,'color','white');
    set(gca,'position',[0.15 0.2 0.7 0.7]);
    set(gca,'Box','on');
    set(gca,'TickDir','out','TickLength',[.02 .02]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    set(gca,'XGrid','off','YGrid','on');
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(gca,'YDir','normal')
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',22);
    set(fh2,'LineWidth',3.0);
    set(fh2,'MarkerSize',20.0);
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
    labX = xlabel(xlab,'FontSize',22,'FontName','Helvetica');
    labY = ylabel(ylab,'FontSize',22,'FontName','Helvetica');
    %set(labX,'interpreter','latex');
    %set(labY,'interpreter','latex');
    
	title(gtitle,'FontSize',22,'FontName','Helvetica');

    if ( strcmp(flag,'eps') )
        saveas(gcf,gname,'epsc2');
        %gname = [gname, '.eps'];
    end


return
% -----------------------------------------------------------------
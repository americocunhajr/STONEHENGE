% -----------------------------------------------------------------
%  myplot_params.m
% ----------------------------------------------------------------- 
%  This function deals with plotting and set-up.
% ----------------------------------------------------------------- 
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
%  programmers: 
%        João Peterson (ligier.peterson@gmail.com)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Oct 20, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = myplot_params(X, Y, X2, Y2, gname, flag)
    fig = figure('NumberTitle','off');
    ax1 = axes('Position',[0.12 0.17 0.7 0.7]);
    plot(ax1, X, Y, 'b')
    xlabel(ax1,'time', 'FontSize', 16, 'FontName', 'Helvetica');
    ylabel(ax1,'displacement', 'FontSize', 16, 'FontName', 'Helvetica');
    set(ax1,'Box','on');
    set(ax1,'TickDir','out','TickLength',[.02 .02]);
    set(ax1,'XMinorTick','on','YMinorTick','on');
    set(ax1,'XGrid','off','YGrid','on');
    set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax1,'FontName','Helvetica');
    set(ax1,'FontSize',14);
    
    ax2 = axes('Position',[0.65 0.65 0.28 0.28]);
    plot(ax2, X2, Y2, 'b')
    set(ax2,'Box','on');
    set(ax2,'TickDir','out','TickLength',[.02 .02]);
    set(ax2,'XMinorTick','on','YMinorTick','on');
    set(ax2,'XGrid','off','YGrid','on');
    set(ax2,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax2,'FontName','Helvetica');
    set(ax2,'FontSize',14);
    
    if (strcmp(flag, 'png'))
        saveas(gcf, gname, 'png')
        
    end
    
end
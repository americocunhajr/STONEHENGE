% -----------------------------------------------------------------
%  plot_3D_bifurcation.m
% ----------------------------------------------------------------- 
%  This function plot the bifurcation diagram by graph 3D.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jp.norenberg@unesp.br)
%
%  last update: Oct 20, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function fig = plot_3D_bifurcation(bifurc_inf,name_file)

    % check number of arguments
     if nargin > 2
         error('Too many inputs.')
     elseif nargin < 2
         error('Too few inputs.')
     end 

    % bifurcation input parameters
     Bifurc_up   = bifurc_inf.Bifurc_up; 
     Bifurc_down = bifurc_inf.Bifurc_down;
     Omega_rang  = bifurc_inf.Omega_rang;
     f_rang_up   = bifurc_inf.f_rang_up;
     f_rang_down = bifurc_inf.f_rang_down;

    % generate 2D mesh grid
     [X1,Y1] = meshgrid(f_rang_up,Omega_rang);
     [X2,Y2] = meshgrid(f_rang_down,Omega_rang);

    % color map
     colorMap = winter(5);
     colorMap2= autumn(5);


    fig = figure(1);
    set(fig, 'Position', [300 50 900 630])
    % loop for plotting each bifurcation
     for k = 1:5
         plot3(X1(k,:),Y1(k,:),Bifurc_up(:,:,k),'.','color',colorMap(k,:), ...
             'MarkerSize',4.2)
         hold on
         plot3(X2(k,:),Y2(k,:),Bifurc_down(:,:,k),'.','color',colorMap2(k,:),...
             'MarkerSize',4.2)
     end

    % set-up plot
     set(gca,'TickDir','out','TickLength',[.02 .02]);
     grid on

     view(40,60)

     xlabel('amplitude','FontSize',15,'FontName','Helvetica');
     ylabel('frequency','FontSize',15,'FontName','Helvetica' );
     zlabel('voltage','FontSize',15,'FontName','Helvetica');

     yh = get(gca,'YLabel'); 
     set(yh, 'Units', 'Normalized')
     pos = get(yh, 'Position');
     set(yh, 'Position',pos.*[1,1,1],'Rotation',35)
     xh = get(gca,'XLabel'); 
     set(xh, 'Units', 'Normalized')
     pos = get(xh, 'Position');
     set(xh, 'Position',pos.*[1,1,1],'Rotation',-25)
 
     xlim([.01 0.3])
     yticks(Omega_rang)
     %ylim([0.64 0.96])
     zlim([-3.5 2.5])

     saveas(fig,name_file,'eps');
end
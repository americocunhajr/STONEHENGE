
% -----------------------------------------------------------------
%  plot_ce_animation_video.m
%
%  This functions plots an animation of cross-entropy method
%  for optimization
%
%  input:
%  x1      - x data vector
%  y1      - y data vector
%  itervec - iteration index vector
%  xref    - reference value for x
%  yref    - reference value for y
%  vtitle  - video title
%  vname   - video file name
%  xlab    - x axis label
%  ylab    - y axis label
%  xmin    - x axis minimum value
%  xmax    - x axis maximum value
%  ymin    - y axis minimum value
%  ymax    - y axis maximum value
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Feb 10, 2016
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function mov = plot_ce_animation_video(x1,y1,itervec,...
                                 xref,yref,vtitle,vname,...
                                 xlab,ylab,xmin,xmax,ymin,ymax)
    
    % check number of arguments
    if nargin < 13
        error('Too few inputs.')
    elseif nargin > 13
        error('Too many inputs.')
    end

    % check arguments
    if size(x1) ~= size(y1)
        error('x1 and y1 must be same dimensions')
    end
    
    % coordinates of legend
    xiter = xmin + 0.1*(xmax-xmin);
    yiter = ymin + 0.9*(ymax-ymin);
    
    % number of iterations
    Niter = length(itervec);
    
    % open video object
    writerObj = VideoWriter(vname);
    writerObj.FrameRate = 1;
    open(writerObj);

    % preallocate memory for movie frames
    nFrames = length(1:Niter);
    mov(1:nFrames) = struct('cdata',[], 'colormap',[]);

    % loop to construct the video
    for n=1:Niter
        
        % initialize video frame
        fig = figure(1000);
        
        % define frame properties
        set(gcf,'color','white');
        set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3],'zColor',[.3 .3 .3]);
        set(gca,'FontName','Helvetica');
        set(gca,'FontSize',16);
        set(gca,'Box','on');
        
        % draw trailer base right arm
        fh0 = plot(x1(:,n),y1(:,n),'xb');
        
        hold on
        
        fh1 = plot(xref,yref,'xr');
        
        set(fh0,'LineWidth',2.0);
        set(fh0,'MarkerSize',10.0);
        set(fh0,'MarkerFaceColor','w');
        set(fh0,'MarkerEdgeColor','b');
        set(fh1,'LineWidth',5.0);
        set(fh1,'MarkerSize',20.0);
        set(fh1,'MarkerFaceColor','w');
        set(fh1,'MarkerEdgeColor','r');
        
        % define labels
        xlabel(xlab,'FontSize',22,'FontName','Helvetica');
        ylabel(ylab,'FontSize',22,'FontName','Helvetica');
        
        % define video frame limits
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        
        % display time counter
        text(xiter,yiter,['iteration = ',num2str(itervec(n))],...
                          'Color','k',...
                          'FontName','Helvetica',...
                          'FontSize',12);

        hold off
        
        % define video title
        title(vtitle,'FontSize',20,'FontName','Helvetica');

        % save movie frame
        mov(n) = getframe(gcf);
        writeVideo(writerObj,mov(n));
        
        pause(1.5)
    end
    
    % close movie object
    close(writerObj);
return
% -----------------------------------------------------------------

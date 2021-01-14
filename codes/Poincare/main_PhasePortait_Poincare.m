% -----------------------------------------------------------------
% main_PhasePortait_Poincare.m
% ----------------------------------------------------------------- 
%  This is the main file for a program which response the phase 
%  portait and Poincare section for the piezo-magneto-elastic beam.
% ----------------------------------------------------------------- 
%  programmers: 
%        João Pedro Norenberg (jpcvalese@gmail.com)
%        Americo Cunha (americo.cunha@uerj.br)
%
%  last update: Dec 20, 2020
% -----------------------------------------------------------------
%% Processing
clc
clear
close all

% program header
% -----------------------------------------------------------
disp(' ---------------------------------------------------')
disp('         Phase Portait and Poincare Section         ')
disp('            Bistable EH (linear coupling)           ')
disp('                                                    ')
disp(' by                                                 ')
disp(' João Pedro Norenberg (Unesp)                       ')
disp(' jpcvalese@gmail.com                                ')
disp(' ---------------------------------------------------')
% -----------------------------------------------------------

% physical parameters:
beta  = 0;         % without nonlinear electromechanical coupling

% range of amplitude excitation
frang = [0.041 0.06 0.083 0.091 0.105 0.115 0.147 0.200 0.25];

% Loop to compute Phase Portait and Poincare Section
for k = 1:length(frang)
    f = frang(k);
    [disp(:,k),velo(:,k)] = phase_portait(f,beta);
    [poincare_disp(:,k),poincare_velo(:,k)] = poincare(f,beta);
end

% Save struct variables
Phase_Poincare.poincare_disp = poincare_disp;
Phase_Poincare.poincare_velo = poincare_velo;
Phase_Poincare.disp          = disp;
Phase_Poincare.velo          = velo;
Phase_Poincare.frang         = frang;
%% Ploting

fh1 = figure(1);
set(fh1, 'Position', [488 342 600 360])

pos1 = [.1 .23 0.8 0.6];
subplot('Position', pos1);

set(gcf,'color','w');
colorMap = [1   0    1 ; 0 0 0; 0    1   0 ; 0 1 1 ;1 0.5 0; ...
            0 0.75 0.62; 0 0 1; 0.5 0.5 0.5; 1 0 0];
                
for k = 1:length(frang)
    plot3(k*ones(length(disp(2/3*end:end,k)),1),disp(2/3*end:end,k), ...
          velo(2/3*end:end,k),'.','color',colorMap(k,:),'MarkerSize',2)
    hold on
    plot3(k*ones(length(poincare_disp(:,k)),1),poincare_disp(:,k), ...
        poincare_velo(:,k),'.','color',[0 0 .5],'MarkerSize',10)
    hold on
end

view(18,50)
xticklabels({'0.041'  , '0.060' , '0.083' , '0.091' , ' 0.105 ' , ...
             '0.115 ' , '0.147' , '0.200' , '0.250' })

xlabel('amplitude of excitation','FontSize',15,'FontName','Helvetica');
ylabel('displacement','FontSize',15,'FontName','Helvetica' );
zlabel('velocity','FontSize',15,'FontName','Helvetica');
yh = get(gca,'YLabel'); 
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,1,1],'Rotation',44)
xh = get(gca,'XLabel'); 
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,1,1],'Rotation',-5)
zlim([-1.5 1.5])
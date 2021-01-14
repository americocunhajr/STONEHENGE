% -----------------------------------------------------------------
%  plot_sobol.m
%
%  This function plots the Sobol Indecies to the piezo-magneto-
%  elastic beam.
% ----------------------------------------------------------------- 
%  programmer: João Pedro C V Norenberg
%              jpcvalese@gmail.com
%
%  last update: Dez 13, 2019
% -----------------------------------------------------------------

function fig = plot_sobol(mySobolAnalysis,order,method)

title_name = ['Sobol Index based on ', method];

Num_par = {'$\xi$';'$\chi$';'$\lambda$';'$\kappa$';'$\Omega$';'$f$'};

color = [0.5 0.5 0.5];

if order == 1
    
    if contains(method,'MC')
        y_mc = mySobolAnalysis.Results.FirstOrder;
        error_mc = mySobolAnalysis.Results.Bootstrap.FirstOrder.ConfLevel;

        y = [y_mc];

    elseif contains(method,'PCE')
        y_pce = mySobolAnalysis.Results.FirstOrder;

        y = [y_pce];
    end

    fig = figure();
    b = bar(y);
    hold on
    b(1).FaceColor = color; 
    b(1).EdgeColor = color;

    if contains(method,'MC')
        er = errorbar(y_mc,error_mc,'LineWidth',2);   
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end

    label_x = Num_par'; 
    label_y = {'First Order'};
    ylim([ 0 1]);
    
elseif order ==2
    [Valu_max, Varldx_max] = maxk(mySobolAnalysis.Results.AllOrders{1,2},5);

    Varldx_numC = mySobolAnalysis.Results.VarIdx{2, 1}(Varldx_max,:);

    for i = 1:5
        Par_plot(i,1) = Num_par(Varldx_numC(i,1));
        Par_plot(i,2) = Num_par(Varldx_numC(i,2));
        Name_par(i,1) = strcat({'$'},Par_plot(i,1),{' '},Par_plot(i,2),{'$'});
    end
    
    y = Valu_max;
    label_x = Name_par';
    label_y = {'Second Order'};
    
    % ------------------------------------------------------------------
    fig = figure();
    b = bar(y);
    hold on
    b(1).FaceColor = color; 
    b(1).EdgeColor = color;
    
    ylim([0 max(y)*1.2])
end
set(gca, 'XTickLabel',label_x, 'XTick',1:numel(label_x),...
        'FontName','Helvetica','FontSize',15,'linewidth',1.2,...
        'TickLabelInterpreter','latex');    

ax = gca; xrule = ax.XAxis; xrule.FontSize = 20; 

set(gca,'position',[0.2 0.2 0.7 0.7]);
set(gca,'Box','on');
set(gca,'TickDir','out','TickLength',[.02 .02]);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'XGrid','on','YGrid','on');

ylabel(label_y,'FontSize',22,'FontName','Helvetica');
title(title_name,'FontSize',15,'FontName','Helvetica');
hold off
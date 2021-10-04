function fig = plot_time_series(X, Y, xlab, ylab, gname)
    Njump1 = 1; % number of steps to jump
    Njump2 = 1;
    Njump3 = 1;

    Ndt = length(X);
    Nss = round(0.9*Ndt);
    Nts1 = round(0.05*Ndt);
    Nts2 = round(0.15*Ndt);
    
    [K, flag] = z1test(Y(Nts1:60:round(0.5*Ndt)));
    
    X2 = X(Nss:Njump2:end);
    Y2 = Y(Nss:Njump2:end);
    
    fig = figure('NumberTitle','off');
    ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
    plot(ax1, X, Y, 'b')
    hold on
    pos1 = X(Nss);
    pos2 = min(Y(Nss:Ndt))-0.1;
    pos3 = (X(end)-X(Nss))*0.2;
    pos4 = max(Y(Nss:Ndt))-min(Y(Nss:Ndt))+0.2;
    rectangle('position',[pos1, pos2, pos3, pos4],...
                  'EdgeColor', 'red', 'LineWidth', 2.0);

    ax2 = axes('Position',[0.65 0.70 0.25 0.21]);
    plot(ax2, X2, Y2, 'b')
    hold on
    xlim([pos1 pos1+pos3]);
    ylim([pos2 pos2+pos4]);
    set(ax2,'XColor', 'red');
    set(ax2,'YColor', 'red');
    set(ax2,'LineWidth', 2.0);
    set(ax2,'Box','on');
    set(ax2,'XTickLabel','');
    set(ax2,'YTickLabel','');
    hold on
    
    if K > 0.8
        X3 = X(Nts1:Njump3:Nts2);
        Y3 = Y(Nts1:Njump3:Nts2);
    
        pos1 = X(Nts1);
        pos2 = min(Y(Nts1:Nts2))-0.1;
        pos3 = (X(Nts2)-X(Nts1))*0.2;
        pos4 = max(Y(Nts1:Nts2))-min(Y(Nts1:Nts2))+0.2;
        rectangle(ax1, 'position',[pos1, pos2, pos3, pos4],...
                      'EdgeColor', 'yellow', 'LineWidth', 2.0);

        ax3 = axes('Position',[0.25 0.70 0.25 0.21]);
        plot(ax3, X3, Y3, 'b')
        hold on
        xlim([pos1 pos1+pos3]);
        ylim([pos2 pos2+pos4]);
        set(ax3,'XColor', 'yellow');
        set(ax3,'YColor', 'yellow');
        set(ax3,'LineWidth', 2.0);
        set(ax3,'Box','on');
        set(ax3,'XTickLabel','');
        set(ax3,'YTickLabel','');
        hold off 
    end
    
    ylim(ax1,[-3 3]);
    set(gcf,'color','white');
    set(ax1,'Box','on');
    set(ax1,'TickDir','out','TickLength',[.02 .02]);
    set(ax1,'XMinorTick','on','YMinorTick','on');
    set(ax1,'XGrid','off','YGrid','on');
    set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
    set(ax1,'FontName','Helvetica');
    set(ax1,'FontSize',14);
    xlabel(ax1,xlab, 'FontSize', 18, 'FontName', 'Helvetica');
    ylabel(ax1,ylab, 'FontSize', 18, 'FontName', 'Helvetica');
    hold off

    saveas(gcf, gname, 'png')
    %close(fig);
end


%% plot displacement
%{
number = 3;

filename = "f_" + number + ".mat";
load(filename);

Njump1 = 1; % number of steps to jump
Njump2 = 1;
Njump3 = 1;

Ndt = length(time);
Nss = round(0.9*Ndt);
Nts1 = round(0.05*Ndt);
Nts2 = round(0.1*Ndt);

Qdisp = y(:,1);

X = time(1:Njump1:end);
Y = Qdisp(1:Njump1:end);
X2 = time(Nss:Njump2:end);
Y2 = Qdisp(Nss:Njump2:end);
X3 = time(Nts1:Njump3:Nts2);
Y3 = Qdisp(Nts1:Njump3:Nts2);

fig1 = figure('NumberTitle','off');
ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
plot(ax1, X, Y, 'b')
hold on
pos1 = X(Nss);
pos2 = min(Y(Nss:Ndt))-0.2;
pos3 = (X(end)-X(Nss))*0.2;
pos4 = abs(max(Y(Nss:Ndt)))+abs(min(Y(Nss:Ndt)))+0.4;
rectangle('position',[pos1, pos2, pos3, pos4],...
              'EdgeColor', 'red', 'LineWidth', 2.0);
          
ax2 = axes('Position',[0.65 0.70 0.25 0.21]);
plot(ax2, X2, Y2, 'b')
hold on
xlim([pos1 pos1+pos3]);
ylim([pos2 pos2+pos4]);
set(ax2,'XColor', 'red');
set(ax2,'YColor', 'red');
set(ax2,'LineWidth', 2.0);
set(ax2,'Box','on');
set(ax2,'XTickLabel','');
set(ax2,'YTickLabel','');
hold on

%{
pos1 = X(Nts1);
pos2 = min(Y(Nts1:Nts2))-0.05;
pos3 = (X(Nts2)-X(Nts1))*0.2;
pos4 = abs(max(Y(Nts1:Nts2)))+abs(min(Y(Nts1:Nts2)))+0.1;
rectangle(ax1, 'position',[pos1, pos2, pos3, pos4],...
              'EdgeColor', 'yellow', 'LineWidth', 2.0);
          
ax3 = axes('Position',[0.25 0.70 0.25 0.21]);
plot(ax3, X3, Y3, 'b')
hold on
xlim([pos1 pos1+pos3]);
ylim([pos2 pos2+pos4]);
set(ax3,'XColor', 'yellow');
set(ax3,'YColor', 'yellow');
set(ax3,'LineWidth', 2.0);
set(ax3,'Box','on');
set(ax3,'XTickLabel','');
set(ax3,'YTickLabel','');
hold off
%}

ylim(ax1,[-3 3]);
set(gcf,'color','white');
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
xlabel(ax1,'time', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel(ax1,'displacement', 'FontSize', 18, 'FontName', 'Helvetica');
hold off

gname = ['disp_f_case_', num2str(number)];
saveas(gcf, gname, 'png')
close(fig1);
%}

%% plot voltage
%{
number = 9;

filename = "f_" + number + ".mat";
load(filename);

Njump1 = 1; % number of steps to jump
Njump2 = 1;
Njump3 = 1;

Ndt = length(time);
Nss = round(0.9*Ndt);
Nts1 = round(0.05*Ndt);
Nts2 = round(0.15*Ndt);

Qvolt = y(:,3);

X = time(1:Njump1:end);
Y = Qvolt(1:Njump1:end);
X2 = time(Nss:Njump2:end);
Y2 = Qvolt(Nss:Njump2:end);
X3 = time(Nts1:Njump3:Nts2);
Y3 = Qvolt(Nts1:Njump3:Nts2);

fig1 = figure('NumberTitle','off');
ax1 = axes('Position',[0.15 0.17 0.7 0.6]);
plot(ax1, X, Y, 'b')
hold on
pos1 = X(Nss);
pos2 = min(Y(Nss:Ndt))-0.1;
pos3 = (X(end)-X(Nss))*0.2;
pos4 = abs(max(Y(Nss:Ndt)))+abs(min(Y(Nss:Ndt)))+0.2;
rectangle('position',[pos1, pos2, pos3, pos4],...
              'EdgeColor', 'red', 'LineWidth', 2.0);
          
ax2 = axes('Position',[0.65 0.70 0.25 0.21]);
plot(ax2, X2, Y2, 'b')
hold on
xlim([pos1 pos1+pos3]);
ylim([pos2 pos2+pos4]);
set(ax2,'XColor', 'red');
set(ax2,'YColor', 'red');
set(ax2,'LineWidth', 2.0);
set(ax2,'Box','on');
set(ax2,'XTickLabel','');
set(ax2,'YTickLabel','');
hold on

%{
pos1 = X(Nts1);
pos2 = min(Y(Nts1:Nts2))-0.1;
pos3 = (X(Nts2)-X(Nts1))*0.2;
pos4 = abs(max(Y(Nts1:Nts2)))+abs(min(Y(Nts1:Nts2)))+0.2;
rectangle(ax1, 'position',[pos1, pos2, pos3, pos4],...
              'EdgeColor', 'yellow', 'LineWidth', 2.0);
          
ax3 = axes('Position',[0.25 0.70 0.25 0.21]);
plot(ax3, X3, Y3, 'b')
hold on
xlim([pos1 pos1+pos3]);
ylim([pos2 pos2+pos4]);
set(ax3,'XColor', 'yellow');
set(ax3,'YColor', 'yellow');
set(ax3,'LineWidth', 2.0);
set(ax3,'Box','on');
set(ax3,'XTickLabel','');
set(ax3,'YTickLabel','');
hold off
%}

ylim(ax1,[-3 3]);
set(gcf,'color','white');
set(ax1,'Box','on');
set(ax1,'TickDir','out','TickLength',[.02 .02]);
set(ax1,'XMinorTick','on','YMinorTick','on');
set(ax1,'XGrid','off','YGrid','on');
set(ax1,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax1,'FontName','Helvetica');
set(ax1,'FontSize',14);
xlabel(ax1,'time', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel(ax1,'voltage', 'FontSize', 18, 'FontName', 'Helvetica');
hold off

gname = ['volt_f_case_', num2str(number)];
saveas(gcf, gname, 'png')
%close(fig1);
%}
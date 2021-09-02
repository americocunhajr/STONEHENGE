% OGY Method

clc
clear all
close all

% Piezomagbeam parameters
ksi    = 0.01;  
chi    = 0.05;  
lambda = 0.05;  
kappa  = 0.5;   
f0     = 0.083; 
Omega  = 0.8;
phys_param = [ksi chi f0 Omega lambda kappa];

% code parameters
N_cycles = 20000;
period = 2*pi/Omega;
N_sampling = 100;
t1 = 0 + N_cycles*period;

time_x0 = 0:period:t1;

tol = 0.05;

% saving memory for x and f vectors
x = zeros(3,N_cycles*N_sampling+1);
f = zeros(N_cycles,1);

% initial conditions
x0 = zeros(3,N_cycles);
x0(1,1) = 1; x0(2,1) = 0; x0(3,1) = 0;
x(:,1) = x0(:,1);

% the values of f will be changed in the control
f(:) = f0;


% flags
orbit_flag = 0;
control_flag = 0;

opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

for n=1:N_cycles
    tspan = (n-1)*period : period/N_sampling : n*period;%-period/N_sampling; 
    [time_aux, y] = ode45(@(t,y)piezomagbeam(t,y,phys_param),tspan,x0(:,n),opt);
    y = y';
    
    % update x0 vector
    x0(:,n+1) = y(:,end);
    
    % update x and time vectors
    x(:,(n-1)*N_sampling+2:n*N_sampling+1) = y(:,2:end);
    time((n-1)*N_sampling+1:n*N_sampling+1) = time_aux(:);
    
    % orbit starting at n-1 = 2-period orbit
    % orbit starting at n   = 1-period orbit
    orbit_start = n;
        
    % define orbit to stabilize the dynamic
    if n > 1 && orbit_flag == 0
        % distance between two points in the Poincaré map
        d = norm(x0(:,orbit_start) - x0(:,n+1));
        
        if d < tol
           % orbit found
           orbit_flag = 1;
            
           % save last period evolution
           if orbit_start == n-1
               % for 2-period orbit it is necessary to save 2 cycles
               orbit(:,:) = [buffer(:,:) y(:,2:end)];
           else
               % single-period orbit
               orbit(:,:) = y(:,2:end);
           end
           
        end
    
    elseif orbit_flag == 1
        % with orbit already selected, check if next point is near any
        % point of the orbit
        
        d_vec = sqrt((x0(1,n+1) - orbit(1,:)).^2 + ...
                     (x0(2,n+1) - orbit(2,:)).^2 + ...
                     (x0(3,n+1) - orbit(3,:)).^2);
        d = min(d_vec);
                
        if d < tol
            % flag indicating in which point the control has been done
            if control_flag == 0
                control_flag = N_sampling*(n+1);


                % find position in d_vec
                pos = find(d_vec == d);

                % if the nearest point is the last of the orbit sequence, seek
                % the first, otherwise, pick the next
                if pos == length(d_vec)
                    pos = 1;
                else
                    pos = pos+1;
                end

                % calculate the impulse
                J = [        0                 1          0     ;
                    0.5*(1.0 - x(1,n+1)^2)   -2*ksi      chi    ;
                             0               -kappa     -lambda ;];

                C = [0 ; cos(Omega*time_x0(n+1));  0 ];

                K(1) = 0.3*sqrt(3)*(1-norm(J))/norm(C);
                K(2) = 0.3*sqrt(3)*(1-norm(J))/norm(C);
                K(3) = 0.3*sqrt(3)*(1-norm(J))/norm(C);

                d1 = orbit(1,pos) - x(1,n+1);
                d2 = orbit(2,pos) - x(2,n+1);
                d3 = orbit(3,pos) - x(3,n+1);

                f(n+1) = -K(1)*d1 -K(2)*d2 -K(3)*d3 + f(n);
            else
                f(n+1) = 0.083;
            end
                
                phys_param = [ksi chi f(n+1) Omega lambda kappa];
        else
            f(n+1) = 0.083;
        end
    
    end
    
    if orbit_start == n-1
        buffer(:,:) = y(:,2:end);
    end
    
end

% Plotting the figures

% voltage time series
gname = ['OGY_volt_f083_[1,0,0]'];
plot_time_series(time,x(3,:),'time','voltage',gname);

% phase space
Ndt = length(x(1,:));
Nss = round(0.9*Ndt);

fig2 = figure('NumberTitle','off');
ax = axes('Position',[0.15 0.17 0.7 0.7]);
%plot(ax,x(1,1:control_flag-100), x(2,1:control_flag-100), '-b')
plot(ax,x(1,1:50000), x(2,1:50000), '-b')
hold on
plot(ax,x(1,Nss:Ndt), x(2,Nss:Ndt), '-r')


set(gcf,'color','white');
set(ax,'Box','on');
set(ax,'TickDir','out','TickLength',[.02 .02]);
set(ax,'XMinorTick','on','YMinorTick','on');
set(ax,'XGrid','off','YGrid','on');
set(ax,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax,'FontName','Helvetica');
set(ax,'FontSize',14);
xlabel(ax,'displacement', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel(ax,'velocity', 'FontSize', 16, 'FontName', 'Helvetica');

% Poincaré map
fig3 = figure('NumberTitle','off');
plot(x0(1,1:control_flag/N_sampling), x0(2,1:control_flag/N_sampling), '.b')
hold on
plot(x0(1,Nss:Ndt), x0(2,Nss:Ndt), '.r')

ax = axes('Position',[0.15 0.17 0.7 0.7]);
set(gcf,'color','white');
set(ax,'Box','on');
set(ax,'TickDir','out','TickLength',[.02 .02]);
set(ax,'XMinorTick','on','YMinorTick','on');
set(ax,'XGrid','off','YGrid','on');
set(ax,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax,'FontName','Helvetica');
set(ax,'FontSize',14);
xlabel(ax,'displacement', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel(ax,'velocity', 'FontSize', 16, 'FontName', 'Helvetica');

% force applied
fig4 = figure('NumberTitle','off');
plot(time_x0,f,'-b')

ax = axes('Position',[0.15 0.17 0.7 0.7]);
set(gcf,'color','white');
set(ax,'Box','on');
set(ax,'TickDir','out','TickLength',[.02 .02]);
set(ax,'XMinorTick','on','YMinorTick','on');
set(ax,'XGrid','off','YGrid','on');
set(ax,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);
set(ax,'FontName','Helvetica');
set(ax,'FontSize',14);
xlabel(ax,'time', 'FontSize', 16, 'FontName', 'Helvetica');
ylabel(ax,'forcing amplitude', 'FontSize', 16, 'FontName', 'Helvetica');
xlim([0 time(end)]);



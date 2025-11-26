function plot_3clusters_9rods_traj
clc; clear all; clf


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
omega = 36*pi;

n_cluster = 3;              % include the one at the center
D    = 4;
rc   = 3*D;     xcenter = rc;
tplt = (pi/3)/omega; 

time = 2*pi/omega;
iter = round(time/tplt); 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
time = tplt*iter; 
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    cd out_data_3clusters_9rods/traj/
        Xpt = load(name_mp); 
        X1 = Xpt(:,1);    X2 = Xpt(:,2);    X3 = Xpt(:,3);
        U1 = Xpt(:,4);    U2 = Xpt(:,5);    U3 = Xpt(:,6);
    cd ..
    
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    figure(1);clf(1)  % xy plane
        plot3(X1-xcenter,X2,X3,'k.','markersize',6); hold on
        plot3([0 rc 2*rc]-xcenter,[0 0 0],[0 0 0],'r.','markersize',24)
        
        axis equal
        xlim([-2*rc-2*rc rc*n_cluster+rc])
        ylim([-35 35])
        set(gca,'fontsize',20)
        xlabel('$x$', 'Interpreter', 'Latex','Fontsize',22)
        ylabel('$y$', 'Interpreter', 'Latex','Fontsize',22)
        zlabel('$z$', 'Interpreter', 'Latex','Fontsize',22)
        box on
        grid on
        view(2)
        pause(0.05)


    
    
    
    
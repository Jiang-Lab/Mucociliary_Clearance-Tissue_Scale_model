function plot_fig15_3cluster_varyD_nrods_mp_1_shift_center
clc; clear all; clf


cd out_data_stokeslet_mp_split
omega = 36*pi;

n_cluster = 3;              % include the one at the center
D    = 4;
rc   = 3*D;     xcenter = rc;
tplt = (pi/3)/omega; 

time = 0; 
% time = 50*pi/omega;
time = 100*pi/omega;
time = 400*pi/omega;


iter = round(time/tplt); 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    time = tplt*iter; 
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    
    cd out_data_1
    Xpt = load(name_mp); 
    X1_1 = Xpt(:,1);    X2_1 = Xpt(:,2);    X3_1 = Xpt(:,3);
    U1_1 = Xpt(:,4);    U2_1 = Xpt(:,5);    U3_1 = Xpt(:,6);
    cd ..
    
    cd out_data_2
    Xpt = load(name_mp); 
    X1_2 = Xpt(:,1);    X2_2 = Xpt(:,2);    X3_2 = Xpt(:,3);
    U1_2 = Xpt(:,4);    U2_2 = Xpt(:,5);    U3_2 = Xpt(:,6);
    cd ..
    
    cd out_data_3
    Xpt = load(name_mp); 
    X1_3 = Xpt(:,1);    X2_3 = Xpt(:,2);    X3_3 = Xpt(:,3);
    U1_3 = Xpt(:,4);    U2_3 = Xpt(:,5);    U3_3 = Xpt(:,6);
    cd ..
    
    cd out_data_4
    Xpt = load(name_mp); 
    X1_4 = Xpt(:,1);    X2_4 = Xpt(:,2);    X3_4 = Xpt(:,3);
    U1_4 = Xpt(:,4);    U2_4 = Xpt(:,5);    U3_4 = Xpt(:,6);
    cd ..
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>>  
    figure(1);clf(1)  % xy plane
        plot3(X1_1-xcenter,X2_1,X3_1,'k.','markersize',6); hold on
        plot3(X1_2-xcenter,X2_2,X3_2,'r.','markersize',6); hold on
        plot3(X1_3-xcenter,X2_3,X3_3,'b.','markersize',6); hold on 
        plot3(X1_4-xcenter,X2_4,X3_4,'g.','markersize',6); hold on 
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

        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    figure(2);clf(2)  % xz plane
        plot3(X1_1-xcenter,X2_1,X3_1,'k.','markersize',6); hold on
        plot3(X1_2-xcenter,X2_2,X3_2,'r.','markersize',6); hold on
        plot3(X1_3-xcenter,X2_3,X3_3,'b.','markersize',6); hold on 
        plot3(X1_4-xcenter,X2_4,X3_4,'g.','markersize',6); hold on 
        plot3([0 rc 2*rc]-xcenter,[0 0 0],[0 0 0],'r.','markersize',24)
        
        axis equal
%         xlim([-2*rc rc*n_cluster+rc])
        xlim([-2*rc-2*rc rc*n_cluster+rc])
        zlim([0 35])
        set(gca,'fontsize',20)
        xlabel('$x$', 'Interpreter', 'Latex','Fontsize',22)
        ylabel('$y$', 'Interpreter', 'Latex','Fontsize',22)
        zlabel('$z$', 'Interpreter', 'Latex','Fontsize',22)
        box on
        grid on
        view(0,0)
        pause(0.05)


   
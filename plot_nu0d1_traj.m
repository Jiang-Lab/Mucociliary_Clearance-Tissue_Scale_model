function plot_nu0d1_traj
clc; clear all; clf


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
omega = 36*pi;
D    = 4;
tplt = (pi/3)/omega; 

time = 2*pi/omega;
iter = round(time/tplt); 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
iter = round(time/tplt); 
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    cd out_data_nu0d1/traj
        Xpt = load(name_mp); 
        X1 = Xpt(:,1);    X2 = Xpt(:,2);    X3 = Xpt(:,3);
        U1 = Xpt(:,4);    U2 = Xpt(:,5);    U3 = Xpt(:,6);
    cd ..

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    figure(1);clf(1)  % xy plane
        plot3(X1,X2,X3,'k.','markersize',6); hold on

        line([-50 -50],[ 50 -50],'color','k')
        line([ 50 -50],[ 50  50],'color','k')
        line([ 50  50],[-50  50],'color','k')
        line([-50  50],[-50 -50],'color','k')
        
        axis equal
        xlim([-100, 100])
        ylim([-100, 100])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; zlabel z; box on
        grid on
        view(2)
        pause(0.05)
    

    
    
    
    
    
    
    
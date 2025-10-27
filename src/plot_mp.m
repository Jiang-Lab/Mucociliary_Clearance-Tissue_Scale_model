function plot_mp
clc; clear all; clf

% mp points
% rod loc + vel

% cd ./out_data


tplt  = pi/8; 
omega = 36*pi;
% time = 0; 
% time = 2*pi/omega; 
% time = 4*pi/omega; 
% time = 8*pi/omega;
time = 500*pi/omega;

% time = 16*pi/omega;
% time = 20*pi/omega;
iter = round(time/tplt); 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    time = tplt*iter; 
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    name_loc_vel  = sprintf('time_loc_vel_%05s',num2str(iter));
    
    Xpt = load(name_mp); 
    X1 = Xpt(:,1);    X2 = Xpt(:,2);    X3 = Xpt(:,3);
    U1 = Xpt(:,4);    U2 = Xpt(:,5);    U3 = Xpt(:,6);
    
    s   = load(name_loc_vel);  
    x  = s(:,1);  y = s(:,2);  z = s(:,3);
    vx = s(:,4); vy = s(:,5); vz = s(:,6);
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>>  
    figure(1);clf(1)
%         plot3(x,y,z,'k.','markersize',25); hold on 
%         quiver3(x,y,z,vx,vy,vz,0.5)  
        plot3(X1,X2,X3,'b.','markersize',1); hold on
%         quiver3(X1,X2,X3,U1,U2,U3); hold on
        
        axis equal
%         xlim([-1.5 1.5])
%         ylim([-1.5 1.5])
% 
%         xlim([-4 4])
%         ylim([-4 4])

%         xlim([-5 5])
%         ylim([-5 5])
        
        xlim([-35 15])
        ylim([-35 15])

        
        zlim([0 20])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; zlabel z; box on
        grid on
%         view(20,40)
%         view(5,5)
        view(2)
%         view(0,0)
        pause(0.05)
    
    pause(0.1)
%         pause

    figure(2);clf(2)
%         plot3(x,y,z,'k.','markersize',25); hold on 
%         quiver3(x,y,z,vx,vy,vz,0.5)  
        plot3(X1,X2,X3,'b.','markersize',1); hold on
%         quiver3(X1,X2,X3,U1,U2,U3); hold on
        
        axis equal
%         xlim([-1.5 1.5])
%         ylim([-1.5 1.5])

        xlim([-4 4])
        ylim([-4 4])
       

        xlim([-5 5])
        ylim([-5 5])
        
        xlim([-3 9])
        zlim([0 1.5])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; zlabel z; box on
        grid on
%         view(20,40)
%         view(5,5)
%         view(2)
        view(0,0)
        pause(0.05)



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
function plot_fig12_1cluster_nrods_mp
clc; clear all; 


cd out_data

omega = 36*pi;
tplt = (pi/3)/omega; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% t=0
    time = 0; iter=0;
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    
    Xpt = load(name_mp); 
    X1 = Xpt(:,1);    X2 = Xpt(:,2);    X3 = Xpt(:,3);
    U1 = Xpt(:,4);    U2 = Xpt(:,5);    U3 = Xpt(:,6);
    n1 = find(X3<1);
    n2 = find(X3>1  & X3<5);
    n3 = find(X3>5  & X3<10);
    
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
time = 0; 
time = 25*pi/omega;
time = 50*pi/omega;
time = 100*pi/omega;
time = 400*pi/omega;

iter = round(time/tplt); 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    time = tplt*iter; 
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    
    Xpt = load(name_mp); 
    X1 = Xpt(:,1);    X2 = Xpt(:,2);    X3 = Xpt(:,3);
    U1 = Xpt(:,4);    U2 = Xpt(:,5);    U3 = Xpt(:,6);
    


    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>>  
    figure(1);clf(1)  % xy plane
        plot3(X1(n1),X2(n1),X3(n1),'k.','markersize',6); hold on
        plot3(X1(n2),X2(n2),X3(n2),'r.','markersize',6); hold on
        plot3(X1(n3),X2(n3),X3(n3),'b.','markersize',6); hold on
        
        axis equal
        xlim([-40 25])
        ylim([-40 25])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; zlabel z; box on
        grid on
        view(2)
        pause(0.05)
    set(gca,'FontSize',20)
    xlabel('$x$', 'Interpreter', 'Latex','Fontsize',22)
    ylabel('$z$', 'Interpreter', 'Latex','Fontsize',22)
        grid on

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    figure(2);clf(2)  % xz plane
        plot3(X1(n1),X2(n1),X3(n1),'k.','markersize',6); hold on
        plot3(X1(n2),X2(n2),X3(n2),'r.','markersize',6); hold on
        plot3(X1(n3),X2(n3),X3(n3),'b.','markersize',6); hold on
        axis equal
        xlim([-40 25])
        zlim([0 35])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; zlabel z; box on
        grid on
        view(0,0)
        pause(0.05)

    set(gca,'FontSize',20)
    xlabel('$x$', 'Interpreter', 'Latex','Fontsize',22)
    ylabel('$z$', 'Interpreter', 'Latex','Fontsize',22)
        grid on


    
    
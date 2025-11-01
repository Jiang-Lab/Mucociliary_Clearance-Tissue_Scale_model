function plot_fig8_groups_in_y_nu
clc; clear all; clf

cd nu0d1_1_uniform

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
omega = 36*pi;
D    = 4;  
tplt = (pi/3)/omega; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% always two groups 
 num_division = 2; 
 z_divide = 10/num_division; 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
time = 0; 
time = 50*pi/omega;
% time = 100*pi/omega;
% time = 200*pi/omega;
time = 400*pi/omega;



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% initial groups at iter=0; 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
    iter = 0; 
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    cd out_data_1
    Xpt = load(name_mp); 
    X1_1 = Xpt(:,1);    X2_1 = Xpt(:,2);    X3_1 = Xpt(:,3);
    cd ..
    
    cd out_data_2
    Xpt = load(name_mp); 
    X1_2 = Xpt(:,1);    X2_2 = Xpt(:,2);    X3_2 = Xpt(:,3);
    cd ..
    
    cd out_data_3
    Xpt = load(name_mp); 
    X1_3 = Xpt(:,1);    X2_3 = Xpt(:,2);    X3_3 = Xpt(:,3);
    cd ..
    
    cd out_data_4
    Xpt = load(name_mp); 
    X1_4 = Xpt(:,1);    X2_4 = Xpt(:,2);    X3_4 = Xpt(:,3);
    cd ..

    X = [X1_1; X1_2; X1_3; X1_4]; 
    Y = [X2_1; X2_2; X2_3; X2_4]; 
    Z = [X3_1; X3_2; X3_3; X3_4]; 

    ng1=[]; ng2=[];
    for k=1:num_division
        ng=[];
        ng = find(Z < z_divide*k+0.1 & Z >= z_divide*(k-1)); 
        if(mod(k,2) == 1)           
            ng1 = [ng1; ng]; 
        else
            ng2 = [ng2; ng];  
        end
    end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
iter = round(time/tplt);
    time = tplt*iter; 
    fprintf('time = %12.5f, iter = %d\n',time,iter)
    
    name_mp       = sprintf('time_mp_%05s',num2str(iter)); 
    
    cd out_data_1
    Xpt = load(name_mp); 
    X1_1 = Xpt(:,1);    X2_1 = Xpt(:,2);    X3_1 = Xpt(:,3);
    cd ..
    
    cd out_data_2
    Xpt = load(name_mp); 
    X1_2 = Xpt(:,1);    X2_2 = Xpt(:,2);    X3_2 = Xpt(:,3);
    cd ..
    
    cd out_data_3
    Xpt = load(name_mp); 
    X1_3 = Xpt(:,1);    X2_3 = Xpt(:,2);    X3_3 = Xpt(:,3);
    cd ..
    
    cd out_data_4
    Xpt = load(name_mp); 
    X1_4 = Xpt(:,1);    X2_4 = Xpt(:,2);    X3_4 = Xpt(:,3);
    cd ..

    X = [X1_1; X1_2; X1_3; X1_4]; 
    Y = [X2_1; X2_2; X2_3; X2_4]; 
    Z = [X3_1; X3_2; X3_3; X3_4]; 

    loc_group(1).X = X(ng1); 
    loc_group(1).Y = Y(ng1); 
    loc_group(1).Z = Z(ng1); 

    loc_group(2).X = X(ng2); 
    loc_group(2).Y = Y(ng2); 
    loc_group(2).Z = Z(ng2); 

    X = []; Y=[]; Z=[]; % release memory
    ng1 = []; ng2 = []; % release memory

    figure(2);clf(2)  % xz plane
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
   
    k = 1; 
    plot3(loc_group(k).X,loc_group(k).Y,loc_group(k).Z,'k.','markersize',1); hold on 
    k = 2; 
    plot3(loc_group(k).X,loc_group(k).Y,loc_group(k).Z,'c.','markersize',1); hold on  
        
        axis equal
        set(gca,'fontsize',30)
        xlabel('$x$', 'Interpreter', 'Latex')
        ylabel('$y$', 'Interpreter', 'Latex')
        zlabel('$z$', 'Interpreter', 'Latex')
        box on
        grid on
        view(0,0)

        xlim([-100 100])
        zlim([0 30])
        pause(0.05)

    
    
    
    
    
    
    
    
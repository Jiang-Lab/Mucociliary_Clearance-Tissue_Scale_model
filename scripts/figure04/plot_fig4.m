function plot_fig4
clc; clear all;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% time-dependent location
% >>>>>>>>>>
omega = 36*pi;
tstart = 0; tend = 500*pi/omega; dt = (2*pi/omega)/50;   % it takes 50 steps for one revolution
time = (3*pi/2)/omega;
iter = round(time/dt); 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% parameter, a single cilium
% >>>>>>>>>>
     L      =   7;               % cilia length \= height
theta1      = pi/3;  
theta2      =  0.8; 
omega       = 36*pi;             % average angular velocity
    Nc      = 10;                % Nc points per cilia
theta_xz    = 0;         


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fluid, cluster
% >>>>>>>>>>
visc    = 1d0; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% cluster + rods + domain Lx, Ly 
% >>>>>>>>>>
Lx = 100; Ly = Lx; 
nu = 1; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% material point,  3D meshgrid, n clusters, n rods
% >>>>>>>>>>
Xpt = []; 

% data_out_1
xrange  = linspace(-50, 50 , 16); 
yrange  = linspace(-50, 50 , 8); 
zrange  = linspace(  0, 20 , 8); 

[xg,yg,zg] = meshgrid(xrange,yrange,zrange); [nx,ny,nz] = size(zg); 

X1 = reshape(xg,nx*ny*nz,1); 
X2 = reshape(yg,nx*ny*nz,1); 
X3 = reshape(zg,nx*ny*nz,1);  

Xpt= [X1; X2; X3]; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% smoothing parameter
% >>>>>>>>>>
e    = 0.1;    
tplt = (pi/3)/omega; wtime=0; wnum = 0;

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % n clusters, n rods, get loc and strength
        name_stokeslet  = sprintf('time_stokeslet_%05s',num2str(iter));
        fprintf('         time = %10.5f, read = %12s \n',time,name_stokeslet);
        rods = load(name_stokeslet); 
        x  = rods(:,1);    y  = rods(:,2);    z  = rods(:,3);
        F1 = rods(:,4);    F2 = rods(:,5);    F3 = rods(:,6);

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    X = [x;y;z]; F = [F1;F2;F3];
    
    pwd
    tic;
        Upt = solve_vel_3D_wallzeq0(Xpt,X,F,visc,e); 
    cpuT = toc;
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % sample points
    % >>>>>>>>>>
    N=length(Xpt)/3;
    X1 = Xpt(1:N);    X2 = Xpt(N+1:2*N);    X3 = Xpt(2*N+1:3*N); 
    U1 = Upt(1:N);    U2 = Upt(N+1:2*N);    U3 = Upt(2*N+1:3*N);
    

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>>  
    figure(2);clf(2)  % xy plane
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
        q = quiver3(X1,X2,X3,U1,U2,U3); hold on
        q.AutoScaleFactor = 3;
        q.Color = 'blue';
        q.LineStyle = '-'; 
        q.AlignVertexCenters = 'on'; 
        plot3(x,y,z,'k.','markersize',10)
        
        axis equal
        xlim([-Lx/2 Lx/2])
        ylim([-Ly/2 Ly/2])

        set(gca,'fontsize',14)
        xlabel('$x$', 'Interpreter', 'Latex','Fontsize',16)
        ylabel('$y$', 'Interpreter', 'Latex','Fontsize',16)
        zlabel('$z$', 'Interpreter', 'Latex','Fontsize',16)
        box on
        grid on
        view(-10,20)
 
    
    
    
    
    
    
    
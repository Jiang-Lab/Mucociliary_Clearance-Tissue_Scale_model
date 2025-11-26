function main_3clusters_9rods_comp_traj
clc; clear all;

addpath(genpath('src'));
mkdir out_data_3clusters_9rods/traj
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% parameter
% >>>>>>>>>>
     L      =   7;               % cilia length \= height
theta1      = pi/3;  
theta2      =  0.8; 
omega       = 36*pi;             % average angular velocity
    Nc      = 10;                % Nc points per cilia
theta_xz    = 0;         

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fluid
% >>>>>>>>>>
visc    = 1d0; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% cluster + rods
% >>>>>>>>>>
D  = 4; 
rc = 3*D;                   % horizontal array
n_cluster = 3;              % include the one at the center
n_cilia_per_cluster = 9; 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
xroot = [0 0 0]; phi0 = 0;       
cs = cos(theta_xz); ss = sin(theta_xz); 

xi = xroot(1) + rc*[0:n_cluster-1]; 
yi = xroot(2)*ones(size(xi)); 
zi = xroot(3)*ones(size(xi)); 

if(n_cilia_per_cluster == 9)
    vecx = [-D  0  D -D  0  D -D  0  D];
    vecy = [-D -D -D  0  0  0  D  D  D];
    vecz = [ 1  1  1  1  1  1  1  1  1];
end

if(n_cilia_per_cluster == 4)
    vecx = [-D  D -D  D];
    vecy = [-D -D  D  D];
    vecz = [ 1  1  1  1];
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% material point,  3D meshgrid, n clusters, n rods
% >>>>>>>>>>
Xpt = []; 

% data_out
xrange  = linspace( -rc, rc*3, 16); 
yrange  = linspace(   -4*D,    4*D, 16); 
zrange  = linspace( 0.05,10 ,16); 


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

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% time-dependent location
% >>>>>>>>>>
tstart = 0; tend = 2*pi/omega; dt = (2*pi/omega)/50;   % it takes 50 steps for one revolution
Niter = (tend-tstart)/dt; Niter = round(Niter); 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for iter=0:Niter+2
    time = iter*dt; 
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % n clusters, n rods, get loc and strength
    cd out_data_3clusters_9rods
        name_stokeslet  = sprintf('time_stokeslet_%05s',num2str(mod(iter,50)));
        fprintf('         time = %10.5f, read = %12s \n',time,name_stokeslet);
        rods = load(name_stokeslet); 
        x  = rods(:,1);    y  = rods(:,2);    z  = rods(:,3);
        F1 = rods(:,4);    F2 = rods(:,5);    F3 = rods(:,6);
    cd ../
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    X = [x;y;z]; F = [F1;F2;F3];

    tic;
        Upt = solve_vel_3D_wallzeq0(Xpt,X,F,visc,e); 
        Xpt = Xpt + dt*Upt; 
    cpuT = toc;
    fprintf('iter = %d/%d, time = %d , cputime = %d \n%d',iter,Niter,time,cpuT)

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % Material points max, min
    % >>>>>>>>>>
    N=length(Xpt)/3;
    X1 = Xpt(1:N);    X2 = Xpt(N+1:2*N);    X3 = Xpt(2*N+1:3*N);
    X3(find(X3 < 0)) = 0; Xpt(2*N+1:3*N)=X3; 
    
    fprintf('     [Xpt_xmax, Xpt_xmin] = [%7.5f  %7.5f]\n',max(max(X1)),min(min(X1)))
    fprintf('     [Xpt_ymax, Xpt_ymin] = [%7.5f  %7.5f]\n',max(max(X2)),min(min(X2)))
    fprintf('     [Xpt_zmax, Xpt_zmin] = [%7.5f  %7.5f]\n',max(max(X3)),min(min(X3)))
    fprintf('     [Xc, Yc, Zc]         = [%7.5f  %7.5f   %7.5f]\n',mean(X1),mean(X2),mean(X3))

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>>  
    if(time>wtime-dt/2)
        N=length(Upt)/3;
        U1 = Upt(1:N);    U2 = Upt(N+1:2*N);    U3 = Upt(2*N+1:3*N);
    
        name  = sprintf('time_mp_%05s',num2str(wnum));
        fprintf('       time = %10.5f, print =%12s \n',time,name);
        cd ./out_data_3clusters_9rods/traj/
        fileID_mp = fopen(name,'w');
        fprintf(fileID_mp,'%12.7e   %12.8e   %12.8e  %12.7e   %12.8e   %12.8e\n',...
            [X1 X2 X3 U1 U2 U3]'); % loop each column
        fclose(fileID_mp);
        cd ../../
        

        wnum  = wnum + 1;
        wtime = wtime + tplt; 
        fprintf('*************************************************************\n')
    end

end


    
    
    
    
    
    
    
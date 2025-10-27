function main_1cluster_1rod
clc; clear all; clf

mkdir out_data
fileID = fopen('./out_data/out_input_parameter','w');
fileID1= fopen('./out_data/out_max_displacement','w');

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% parameter, a single cilium
% >>>>>>>>>>
     L      =   7;               % cilia length \= height
theta1      = pi/3;  
theta2      =  0.8; 
omega       = 36*pi;             % average angular velocity
    Nc      = 10;                % Nc points per cilia
% theta_xz    = -pi/12;            % rotation in xz plane
theta_xz    = 0;  
% [cos -sin  0
%  sin  cos  0
%   0    0   1]

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fluid, cluster
% >>>>>>>>>>
visc    = 1d0; 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% rod roots
% >>>>>>>>>>
xroot = [0 0 0]; phi0 = 0;       
cs = cos(theta_xz); ss = sin(theta_xz); 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% material point,  3D meshgrid, 1 cluster, 1 rod
% >>>>>>>>>>
Xpt = []; 
xrange  = linspace(-15,15,16); 
yrange  = linspace(-15,15,16); 
zrange  = linspace( 0.05,10,16); 

[xg,yg,zg] = meshgrid(xrange,yrange,zrange); [nx,ny,nz] = size(zg); 

X1 = reshape(xg,nx*ny*nz,1); 
X2 = reshape(yg,nx*ny*nz,1); 
X3 = reshape(zg,nx*ny*nz,1);  

Xpt= [X1; X2; X3]; 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% smoothing parameter
% >>>>>>>>>>
e    = 0.1;    
tplt = pi/8; wtime=0; wnum = 0;


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% time-dependent location
% >>>>>>>>>>
tstart = 0; tend = 1000*pi/omega;  dt = (2*pi/omega)/50;   % it takes 50 steps for one revolution 
Niter = (tend-tstart)/dt; Niter = round(Niter); 


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fprintf(fileID,'********************\n'); 
fprintf(fileID,'1 cluster, 1 rod\n');
fprintf(fileID,'********************\n');
fprintf(fileID,'\n');

fprintf(fileID,'L        = %d\n',L);
fprintf(fileID,'theta1   = %d\n',theta1);      
fprintf(fileID,'theta2   = %d\n',theta2); 
fprintf(fileID,'omega    = %d\n', omega);
fprintf(fileID,'Nc       = %d\n', Nc);
fprintf(fileID,'theta_xz = %d\n',theta_xz);
fprintf(fileID,'\n');

fprintf(fileID,'visc     = %d\n',visc);
fprintf(fileID,'\n');

fprintf(fileID,'e        = %d\n',e);
fprintf(fileID,'tplt     = %d, wtime = %d, wnum = %d\n', tplt,wtime,wnum); 
fprintf(fileID,'tstart   = %d, tend = %d,  dt = %d\n', tstart,tend,dt);
fprintf(fileID,'***********************************************************\n');
fprintf(fileID,'\n');


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for iter=0:Niter
% for iter=0:10
% for iter=1:1
    time = iter*dt; 
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % 1 cluster, 1 rod, get loc and vel
    arg   = omega*time; 
      w   = abs(sin(arg/2));
    theta = w*theta1+(1-w)*theta2;

    Lc=L; 
    if( mod(arg,2*pi) <= 2*pi && mod(arg,2*pi) >= pi)
       Lc= length_rod(mod(arg,2*pi),L);
    end
    st=linspace(0,Lc,Nc);  s = st(2:end)';  % remove point on xy-plane
    x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[];
    
            % >>>>>>>>>>>>>>>>>>>>
            % location
            % >>>>>>>>>>>>>>>>>>>>
            xc=s*cos(theta)*cos(phi0);
            yc=s*cos(theta)*sin(phi0);
            zc=s*sin(theta); 

            xt = cos(arg).*xc - sin(arg).*yc; 
            yt = sin(arg).*xc + cos(arg).*yc; 
            zt = zc;
            
            x =  cs * xt - ss * zt + xroot(1);
            y =  yt                + xroot(2);
            z =  ss * xt + cs * zt + xroot(3);
            
            % >>>>>>>>>>>>>>>>>>>>
            % velocity
            % >>>>>>>>>>>>>>>>>>>>
            iarg = mod(arg,2*pi);
            dthetadt= 0.5*omega*cos(iarg/2)*(theta1-theta2); 
            
            dxcdt = -s*sin(theta)*cos(phi0)*dthetadt;
            dycdt = -s*sin(theta)*sin(phi0)*dthetadt; 
            dzcdt =  s*cos(theta)*dthetadt; 

            dxtdt = cos(arg)*dxcdt - sin(arg)*dycdt - omega*sin(arg)*xc - omega*cos(arg)*yc;
            dytdt = sin(arg)*dxcdt + cos(arg)*dycdt + omega*cos(arg)*xc - omega*sin(arg)*yc;
            dztdt = dzcdt; 
            
            vx = cs*dxtdt - ss*dztdt; 
            vy = dytdt;
            vz = ss*dxtdt + cs*dztdt; 
            
%             size(x)
%             size(vx)
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    X = [x;y;z]; U = [vx;vy;vz];
%     size(X)
%     size(U)
    
    tic;
        F   = stokes3d_wallzeq0(X,U,visc,e); 
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
    fprintf(fileID1,'%10.5f  %7.5f  %7.5f  %7.5f  %7.5f  %7.5f  %7.5f \n',...
time,max(max(X1)),min(min(X1)),max(max(X2)),min(min(X2)),max(max(X3)),min(min(X3)));

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>>  
    figure(1);clf(1)
        plot3(x,y,z,'k.','markersize',25); hold on 
%         quiver3(x,y,z,vx,vy,vz,0.5)  
        plot3(X1,X2,X3,'b.','markersize',1); hold on

        axis equal
        xlim([-15 15])
        ylim([-15 15])
        zlim([0 10])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; zlabel z; box on
        grid on
        view(20,40)
%         view(5,5)
%         view(2)
        pause(0.05)

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % write out data
    % >>>>>>>>>> 
    if(time>wtime-dt/2)
        N=length(Upt)/3;
        U1 = Upt(1:N);    U2 = Upt(N+1:2*N);    U3 = Upt(2*N+1:3*N);
    
        name  = sprintf('time_mp_%05s',num2str(wnum));
          fprintf(fileID,'  time = %10.5f, print =%12s \n',time,name);
          fprintf('       time = %10.5f, print =%12s \n',time,name);
        cd ./out_data/
        fileID_mp = fopen(name,'w');
        fprintf(fileID_mp,'%12.7e   %12.8e   %12.8e  %12.7e   %12.8e   %12.8e\n',...
            [X1 X2 X3 U1 U2 U3]'); % loop each column
        fclose(fileID_mp);
        cd ../ 

        
        name  = sprintf('time_loc_vel_%05s',num2str(wnum));
          fprintf(fileID,'  time = %10.5f, print =%12s \n',time,name);
          fprintf('       time = %10.5f, print =%12s \n',time,name);
        cd ./out_data/
        fileID_loc_vel = fopen(name,'w');
        fprintf(fileID_loc_vel,'%12.7e   %12.8e   %12.8e  %12.7e   %12.8e   %12.8e\n',...
            [x y z vx vy vz]'); % loop each column
        fclose(fileID_loc_vel);
        cd ../ 
        

        wnum  = wnum + 1;
        wtime = wtime + tplt; 
        fprintf('*************************************************************\n')
    end
    
    pause(0.1)
%         pause


end


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fclose(fileID);
fclose(fileID1);
return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
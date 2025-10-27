function plot_1cluster_1rod_grid_vel
clc; clear all; clf

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
tstart = 0; 

    time = 0/omega;
%     time = (pi/2)/omega;
%     time = pi/omega;
%     time = (3*pi/2)/omega;
    fprintf('time = %12.5f\n',time)
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % 1 cluster, 1 rod, get loc and vel
    arg   = omega*time; 
      w   = abs(sin(arg/2));
    theta = w*theta1+(1-w)*theta2;
    Lc=L; 
    if( mod(arg,2*pi) <= 2*pi && mod(arg,2*pi) >= pi)
       Lc= length_rod(mod(arg,2*pi),Lc);
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
%        return     
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    X = [x;y;z]; U = [vx;vy;vz];
%     z
    tic;
        F   = stokes3d_wallzeq0(X,U,visc,e); 
        Upt = solve_vel_3D_wallzeq0(Xpt,X,F,visc,e); 
    cpuT = toc;


    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % sample points
    % >>>>>>>>>>
    N=length(Xpt)/3;
    X1 = Xpt(1:N);    X2 = Xpt(N+1:2*N);    X3 = Xpt(2*N+1:3*N);
    X3(find(X3 < 0)) = 0; Xpt(2*N+1:3*N)=X3; 

    U1 = Upt(1:N);    U2 = Upt(N+1:2*N);    U3 = Upt(2*N+1:3*N);
    
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % for visualization
    % >>>>>>>>>> 
    figure(1);clf(1) %     xy-plane 
%         plot3(x,y,z,'k.','markersize',25); hold on 
%         quiver3(x,y,z,vx,vy,vz,0.5)  
%         plot3(X1,X2,X3,'b.','markersize',1); hold on
        plot(x,y,'k.','markersize',25); hold on 
        quiver(X1,X2,U1,U2,10); hold on  
        
        axis equal
        xlim([-25 25])
        ylim([-25 25])
        zlim([0 10])
        set(gca,'fontsize',20)
        xlabel x; ylabel y; box on
        grid on

        
     
    figure(2);clf(2) % xz-plane
%         plot3(x,y,z,'k.','markersize',25); hold on 
%         quiver3(x,y,z,vx,vy,vz,0.5)  
%         plot3(X1,X2,X3,'b.','markersize',1); hold on

        plot(x,z,'k.','markersize',25); hold on 
        quiver(X1,X3,U1,U3,10); hold on  
        
        
        axis equal
        xlim([-25 25])
        ylim([0 15])
        set(gca,'fontsize',20)
        xlabel x; ylabel z; box on
        grid on
        
        



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
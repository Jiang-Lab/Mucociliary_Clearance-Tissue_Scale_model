function plot_fig11right_1cluster_nrods_grid_vel_xzplane
clc; clear all; clf

figure(1);clf(1) %     xz-plane  
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
% rod roots
% >>>>>>>>>>
xroot = [0 0 0]; phi0 = 0;       
cs = cos(theta_xz); ss = sin(theta_xz); 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% parameter, cluster
% >>>>>>>>>>
D  = 4; 
n_cilia_per_cluster = 9;

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
% material point,  3D meshgrid, 1 cluster, 1 rod
% >>>>>>>>>>
Xpt = []; 
xrange  = linspace(-15,15,16); 
yrange  = 0;
zrange  = linspace( 0.05,15,16); 


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
    time = (1*pi/2)/omega;
    time = (pi)/omega;
    time = (3*pi/2)/omega;

    fprintf('time = %12.5f\n',time)
    
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
    
    xx = vecz*xroot(1) + vecx;  
    yy = vecz*xroot(2) + vecy;
    zz = vecz*xroot(3);
    
        x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[];
        for i = 1:n_cilia_per_cluster     % 3x3 location
            
            % >>>>>>>>>>>>>>>>>>>>
            % location
            % >>>>>>>>>>>>>>>>>>>>
            xc=s*cos(theta)*cos(phi0);
            yc=s*cos(theta)*sin(phi0);
            zc=s*sin(theta);   

            xt = cos(arg).*xc - sin(arg).*yc; 
            yt = sin(arg).*xc + cos(arg).*yc; 
            zt = zc;

            % first rotation in xy plane, then shift the root to the
            % desired location
            xtemp =  cs * xt - ss * zt + xx(i);
            ytemp =  yt                + yy(i);
            ztemp =  ss * xt + cs * zt + zz(i);
            
%             plot(xtemp(end),ztemp(end),'ro','markersize',6,'markerfacecolor','r'); hold on 
            
            x = [x; xtemp];    
            y = [y; ytemp]; 
            z = [z; ztemp]; 

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
            
            dxdt = cs*dxtdt - ss*dztdt; 
            dydt = dytdt;
            dzdt = ss*dxtdt + cs*dztdt; 
            
            vx = [vx; dxdt];    
            vy = [vy; dydt]; 
            vz = [vz; dzdt];  
            
        end    
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
    figure(1);clf(1)
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

        quiver(X1,X3,U1,U3,4); hold on  
        plot3(0,0,0,'ro','markersize',8); hold on 
        axis equal
        xlim([-20 20])
        ylim([0 20])
        set(gca,'FontSize',26)
        xlabel('$x$', 'Interpreter', 'Latex','Fontsize',24)
        ylabel('$z$', 'Interpreter', 'Latex','Fontsize',24)
        grid on

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
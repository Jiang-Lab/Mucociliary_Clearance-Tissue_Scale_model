function main_nu0d1_store_strength
clc; clear all; 

addpath(genpath('src'));

mkdir out_data_nu0d1
fileID = fopen('./out_data_nu0d1/out_input_parameter','w');

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
n_cilia_per_cluster = 9; 
D  = 4;  

read_roots = load('input_roots_nu0d1_1');

nu        = read_roots(1,1);
n_cluster = read_roots(1,2); 
xi  = read_roots([2:end],1); 
yi  = read_roots([2:end],2); 
zi =  0*ones(size(xi));

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% angular velocity of each cluster
% >>>>>>>>>>
omega_vec = omega*ones(1,n_cluster);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% phase shift each row 
% >>>>>>>>>>
phi0_shift = 0;
phi0_vec = [1 1 1   % phi vary in x dir
            2 2 2
            0 0 0]*phi0_shift;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cs = cos(theta_xz); ss = sin(theta_xz); 
if(n_cilia_per_cluster == 9)
    vecx = [-D  0  D -D  0  D -D  0  D];
    vecy = [-D -D -D  0  0  0  D  D  D];
    vecz = [ 1  1  1  1  1  1  1  1  1];
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% smoothing parameter
% >>>>>>>>>>
e    = 0.1;    
tplt = (pi/3)/omega; wtime=0; wnum = 0;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% time-dependent location
% >>>>>>>>>>
tstart = 0; tend = 2*pi/omega;  dt = (2*pi/omega)/50;   % it takes 50 steps for one revolution 
Niter = (tend-tstart)/dt; Niter = round(Niter); 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fprintf(fileID,'********************\n'); 
fprintf(fileID,'n clusters random, n rod\n');
fprintf(fileID,'********************\n');
fprintf(fileID,'\n');

fprintf(fileID,'L        = %d\n',L);
fprintf(fileID,'theta1   = %d\n',theta1);      
fprintf(fileID,'theta2   = %d\n',theta2); 
fprintf(fileID,'omega    = %d\n', omega);
for k=1:n_cluster
    fprintf(fileID,'k, omega    = %d  %d\n', k, omega_vec(k));
end
fprintf(fileID,'Nc       = %d\n', Nc);
fprintf(fileID,'theta_xz = %d\n',theta_xz);
fprintf(fileID,'\n');

fprintf(fileID,'visc     = %d\n',visc);
fprintf(fileID,'\n');

fprintf(fileID,'phi0_shift     = %d\n',phi0_shift);
fprintf(fileID,'\n');

fprintf(fileID,'Domain, Lx = %d, Ly = %d\n',Lx,Ly);
fprintf(fileID,'nu (density)  = %5.2f\n',nu);
fprintf(fileID,'D     = %d\n',D);
fprintf(fileID,'n_cluster               = %d\n',n_cluster);
fprintf(fileID,'n_cilia_per_cluster     = %d\n',n_cilia_per_cluster);
fprintf(fileID,'\n');

fprintf(fileID,'e        = %d\n',e);
fprintf(fileID,'tplt     = %d, wtime = %d, wnum = %d\n', tplt,wtime,wnum); 
fprintf(fileID,'tstart   = %d, tend = %d,  dt = %d\n', tstart,tend,dt);
fprintf(fileID,'***********************************************************\n');
fprintf(fileID,'\n');

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for iter=0:Niter+2
    time = iter*dt; 

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % n clusters, n rods, get loc and vel
    x=[]; y=[]; z=[]; vx=[]; vy=[]; vz=[];
    for k=1:n_cluster      % cluster: all cilia roots
    iomega = omega_vec(k); % each cluster has its own velocity

    arg   = iomega*time; 
      w   = abs(sin(arg/2));
    theta = w*theta1+(1-w)*theta2;  
    
    Lc=L; 
    if( mod(arg,2*pi) <= 2*pi && mod(arg,2*pi) >= pi)
       Lc= length_rod(mod(arg,2*pi),L);
    end
%     fprintf('arg = %7.5f, theta = %7.5f, Lc = %7.5f\n',mod(arg,2*pi),theta, Lc)
    st=linspace(0,Lc,Nc);  s = st(2:end)';  % remove point on xy-plane
    
    xx = vecz*xi(k) + vecx;  
    yy = vecz*yi(k) + vecy;
    zz = vecz*zi(k);
    
        for i = 1:n_cilia_per_cluster     % 3x3 location
            % >>>>>>>>>>>>>>>>>>>>
            % location
            % >>>>>>>>>>>>>>>>>>>>
            iphi0 = phi0_vec(i); 

            xc=s*cos(theta)*cos(iphi0);
            yc=s*cos(theta)*sin(iphi0);
            zc=s*sin(theta);   

            xt = cos(arg).*xc - sin(arg).*yc; 
            yt = sin(arg).*xc + cos(arg).*yc; 
            zt = zc;

            % first rotation in xy plane, then shift the root to the
            % desired location
            xtemp =  cs * xt - ss * zt + xx(i);
            ytemp =  yt                + yy(i);
            ztemp =  ss * xt + cs * zt + zz(i);

            x = [x; xtemp];    
            y = [y; ytemp]; 
            z = [z; ztemp]; 
            
            
            % >>>>>>>>>>>>>>>>>>>>
            % velocity
            % >>>>>>>>>>>>>>>>>>>>
            iarg = mod(arg,2*pi);
            dthetadt= 0.5*iomega*cos(iarg/2)*(theta1-theta2); 
            
            dxcdt = -s*sin(theta)*cos(iphi0)*dthetadt;
            dycdt = -s*sin(theta)*sin(iphi0)*dthetadt; 
            dzcdt =  s*cos(theta)*dthetadt; 

            dxtdt = cos(arg)*dxcdt - sin(arg)*dycdt - iomega*sin(arg)*xc - iomega*cos(arg)*yc;
            dytdt = sin(arg)*dxcdt + cos(arg)*dycdt + iomega*cos(arg)*xc - iomega*sin(arg)*yc;
            dztdt = dzcdt; 
            
            dxdt = cs*dxtdt - ss*dztdt; 
            dydt = dytdt;
            dzdt = ss*dxtdt + cs*dztdt; 
            
            vx = [vx; dxdt];    
            vy = [vy; dydt]; 
            vz = [vz; dzdt];  
            
        end
    end
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    X = [x;y;z]; U = [vx;vy;vz];   
    tic;
        F   = stokes3d_wallzeq0(X,U,visc,e); 
    cpuT = toc;
    fprintf('iter = %d/%d, time = %d , cputime = %d \n%d',iter,Niter,time,cpuT)

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % store stokeslet strength
    % >>>>>>>>>> 
    N=length(X)/3;
    F1 = F(1:N);    F2 = F(N+1:2*N);    F3 = F(2*N+1:3*N);
 
        name  = sprintf('time_stokeslet_%05s',num2str(iter));
          fprintf(fileID,'  time = %10.5f, print =%12s \n',time,name);
          fprintf('         time = %10.5f, print =%12s \n',time,name);
        cd ./out_data_nu0d1/
        fileID_strength = fopen(name,'w');
        fprintf(fileID_strength,'%12.7e   %12.8e   %12.8e  %12.7e   %12.8e   %12.8e\n',...
            [x y z  F1 F2 F3]'); % loop each column
        fclose(fileID_strength);
        cd ../ 
    
    pause(0.05)

end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fclose(fileID);

return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

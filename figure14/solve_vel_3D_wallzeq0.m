function Ugrid = solve_vel_3D_wallzeq0(Xgrid,X,F,visc,e)

% This function solves the velocities of specific locations from the 
% strength of the 3D Stokeslet. (Cortez et al (2005))
%
% input variable X saves locations of the Stokeslets. It is a 3N by 1
% vector ([x;y;z]);
% input varibale F saves the strength of the Stokeslet
% input variable e is the regularization parameter.
% input variable Xp saves the locations of interest (the grid)
%
% output variable ugrid saves the computated velocities on the grid. It is
% a 3Np by 1 vector.
%
% Hanliang Guo
% Nov 27 2016

% No slip wall at z=0, modified 2020 June by Ling Xu 


N=length(X)/3;
Ngrid = length(Xgrid)/3;
ugrid = zeros(Ngrid,N);
vgrid = zeros(Ngrid,N);
wgrid = zeros(Ngrid,N);

fx = F(1:N);
fy = F(N+1:2*N);
fz = F(2*N+1:3*N);

xgrid = Xgrid(1:Ngrid);
ygrid = Xgrid(Ngrid+1:2*Ngrid);
zgrid = Xgrid(2*Ngrid+1:end);

X1 = X(1:N);
X2 = X(N+1:2*N);
X3 = X(2*N+1:3*N);

e2= e^2; 
for k = 1:N
        hk = abs(X3(k)); % wall at z=0
%         hk = 0; 
        diffx  = xgrid-X1(k); diffy  = ygrid-X2(k); diffz  = zgrid-X3(k); 
        r2     = diffx.^2  + diffy.^2  + diffz.^2;

        diffxs = xgrid-X1(k); diffys = ygrid-X2(k); diffzs = zgrid+X3(k); 
        r2s    = diffxs.^2 + diffys.^2 + diffzs.^2;
        
        % ***************          
        S11 = func_H1(r2,e2)            - func_H1(r2s,e2) ...
              + func_H2(r2,e2).*diffx.*diffx - func_H2(r2s,e2).*diffxs.*diffxs ...
              + hk^2*(   func_D1(r2s,e2) + func_D2(r2s,e2).*diffxs.*diffxs ) ...
              + 2*hk*(  -func_H2(r2s,e2) .*diffzs ...
                        +0.5*func_D2(r2s,e2).*(-diffxs).*diffzs.*diffxs) ...
              + 2*hk*func_H3(r2s,e2).*(-diffzs); 
          
        S12 = func_H2(r2,e2).*diffy.*diffx - func_H2(r2s,e2).*diffys.*diffxs ...
              + hk^2*(                   + func_D2(r2s,e2) .*diffys.*diffxs ) ...
              + 2*hk*(  +0.5*func_D2(r2s,e2).*(-diffys).*diffzs.*diffxs); 
        
        
        S13 = func_H2(r2,e2).*diffz.*diffx - func_H2(r2s,e2).*diffzs.*diffxs ...
              + hk^2*(                   + func_D2(r2s,e2) .*(-diffzs).*diffxs ) ...
              + 2*hk*(   func_H2(r2s,e2) .*diffxs ...
                        +0.5*func_D2(r2s,e2).*diffzs.*diffzs.*diffxs); 
        
        % ***************
        S21 = func_H2(r2,e2).*diffx.*diffy - func_H2(r2s,e2).*diffxs.*diffys ...
              + hk^2*(                   + func_D2(r2s,e2) .*diffxs.*diffys ) ...
              + 2*hk*(  +0.5*func_D2(r2s,e2).*(-diffxs).*diffzs.*diffys); 
        
        
        S22 = func_H1(r2,e2)            - func_H1(r2s,e2) ...
              + func_H2(r2,e2).*diffy.*diffy - func_H2(r2s,e2).*diffys.*diffys ...
              + hk^2*(   func_D1(r2s,e2) + func_D2(r2s,e2).*diffys.*diffys ) ...
              + 2*hk*(  -func_H2(r2s,e2) .*diffzs ...
                        +0.5*func_D2(r2s,e2).*(-diffys).*diffzs.*diffys) ...
              + 2*hk*func_H3(r2s,e2).*(-diffzs); 
        
        
        S23 = func_H2(r2,e2).*diffz.*diffy - func_H2(r2s,e2).*diffzs.*diffys ...
              + hk^2*(                   + func_D2(r2s,e2) .*(-diffzs).*diffys ) ...
              + 2*hk*(   func_H2(r2s,e2) .*diffys ...
                        +0.5*func_D2(r2s,e2).*diffzs.*diffzs.*diffys); 
        
        % ***************
        S31 = func_H2(r2,e2).*diffx.*diffz - func_H2(r2s,e2).*diffxs.*diffzs ...
              + hk^2*(                   + func_D2(r2s,e2) .*diffxs.*diffzs ) ...
              + 2*hk*(  +0.5*func_D2(r2s,e2).*(-diffxs).*diffzs.*diffzs ...
                        +(func_H3(r2s,e2) - func_H2(r2s,e2)).*(-diffxs) )...
              + 2*hk*func_H3(r2s,e2).*diffxs; 
        
        
        S32 = func_H2(r2,e2).*diffy.*diffz - func_H2(r2s,e2).*diffys.*diffzs ...
              + hk^2*(                   + func_D2(r2s,e2) .*diffys.*diffzs ) ...
              + 2*hk*(  +0.5*func_D2(r2s,e2).*(-diffys).*diffzs.*diffzs ...
                        +(func_H3(r2s,e2) - func_H2(r2s,e2)).*(-diffys) )...
              + 2*hk*func_H3(r2s,e2).*diffys; 
        
        
        S33 = func_H1(r2,e2)            - func_H1(r2s,e2) ...
              + func_H2(r2,e2).*diffz.*diffz - func_H2(r2s,e2).*diffzs.*diffzs ...
              + hk^2*(  -func_D1(r2s,e2) - func_D2(r2s,e2).*diffzs.*diffzs ) ...
              + 2*hk*(  2*func_H2(r2s,e2) .*diffzs ...
                        +0.5*func_D2(r2s,e2).*diffzs.*diffzs.*diffzs ...
                        +(func_H3(r2s,e2) - func_H2(r2s,e2)).*diffzs   );
        
        
        ugrid(:,k) = ugrid(:,k) + (S11*fx(k)+S12*fy(k)+S13*fz(k));
        vgrid(:,k) = vgrid(:,k) + (S21*fx(k)+S22*fy(k)+S23*fz(k));
        wgrid(:,k) = wgrid(:,k) + (S31*fx(k)+S32*fy(k)+S33*fz(k));
end

usurf = sum(ugrid,2);
vsurf = sum(vgrid,2);
wsurf = sum(wgrid,2);

Ugrid = [usurf;vsurf;wsurf]/visc;







function generate_random_number_2D

clc; clear all; 

Lx = 100; 
Ly = 100;  
D=4; 
garea = Lx*Ly; 
larea = 10*10; 
nu    = 0.1; 
num = round(garea*nu/larea);
% num = 10; 
Dmin = 3*D; 


fprintf('nu = %5.2f, num = %d\n', nu,num)
% return
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
% adjust the location by hand
% generate once, then adjust
% >>>>>>>>>>>>>>>>>>>>>>>>>>
   
dum = -Lx/2 + Lx*rand(num,2);
Dist = pdist(dum); 
Z    = squareform(Dist);
[row,col]  = find(Z < Dmin);

count = 0;  
while min(Dist) < Dmin

    count = count + 1;
    
    figure(1);clf(1)
    for k=1:length(row)
        irow = row(k);
        icol = col(k);
        if(icol > irow) % upper triangle matrix
            plot(dum(irow,1),dum(irow,2),'bo','markersize',20); hold on
            plot(dum(icol,1),dum(icol,2),'go','markersize',20); hold on     
        end
    end
            
    xi = dum(:,1);  yi = dum(:,2); theta = linspace(0,2*pi,100); 
    for k = 1:length(xi)
        xx = xi(k) + Dmin*cos(theta);
        yy = yi(k) + Dmin*sin(theta);
        plot(xx,yy,'k-') 
    end
    
    plot(xi,yi,'r.','markersize',20); hold on
    set(gca,'fontsize',20)
    axis equal
    axis([-Lx Lx -Ly Ly])
    xlabel x; ylabel y
    title('before')
%     pause
    pause(0.1)
%     dum

    for k=1:length(row)
        irow = row(k);
        icol = col(k);
        if(icol > irow) % upper triangle matrix
            angle = 2*pi*randn(1,1);  % [0,2pi] 
            dum(icol,1) = dum(irow,1) + Dmin*cos(angle);   
            dum(icol,2) = dum(irow,2) + Dmin*sin(angle);  
        end
    end
%     dum
    
    figure(2);clf(2)
    xi = dum(:,1);
    yi = dum(:,2);
    theta = linspace(0,2*pi,100); 
    for k = 1:length(xi)
        xx = xi(k) + Dmin*cos(theta);
        yy = yi(k) + Dmin*sin(theta);
        plot(xx,yy,'k-'); hold on
    end
    
    plot(xi,yi,'r.','markersize',20); hold on
    set(gca,'fontsize',20)
    axis equal
    axis([-Lx Lx -Ly Ly])
    xlabel x; ylabel y
    title('after')
%     pause
    pause(0.1)

    Dist = pdist(dum); 
    Z    = squareform(Dist);
    [row,col]  = find(Z < Dmin);

    if(count > 200)
       fprintf('cluster num = %d, density = %5.2f, count = %d\n',num,nu,count);
        disp('generating locations, count >=100') 
       return
    end
    
end

        fprintf('generate input_roots,nu = %5.2f, num = %d \n',nu,num)
        dum
        
        name  = sprintf('input_roots');
        fileID_input_roots = fopen(name,'w');
        fprintf(fileID_input_roots,'%12.7f   %d\n',...
            nu,num); % density, number    
        fprintf(fileID_input_roots,'%12.7f   %12.8f\n',...
            dum');   % loop each column
        fclose(fileID_input_roots);


        
        
        
        
        
        
        
        
        
        
        
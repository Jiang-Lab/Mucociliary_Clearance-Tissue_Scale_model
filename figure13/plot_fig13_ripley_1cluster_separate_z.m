function plot_fig13_ripley_1cluster_separate_z

clc
cd ripley_uniform_1cluster_separate_z_xydirs

time_vec  = [1500];
% Area
D  = 4; 
A = (8*D)^2; 

step = 2;
figure(1);clf(1); 
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_hight = 1;  
time_iter=time_vec;
name = sprintf('%s_time%05s_z%05s','data_1cluster',num2str(time_iter),num2str(z_hight)); 
s  = load(name);
rg = s(1:step:end,2); ripley = s(1:step:end,3); 
ripley=smooth(ripley,2);ripley=ripley*A; 

plot(rg,ripley,'ko-','linewidth',2,'markersize',10); hold on
% loglog(rg,ripley,'ko-','linewidth',2,'markersize',10); hold on


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_hight = 3;  
time_iter=time_vec;
name = sprintf('%s_time%05s_z%05s','data_1cluster',num2str(time_iter),num2str(z_hight)); 
s  = load(name);
rg = s(1:step:end,2); ripley = s(1:step:end,3); 
ripley=smooth(ripley,2);ripley=ripley*A; 
loglog(rg,ripley,'r>-','linewidth',2,'markersize',10); hold on


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_hight = 5;  
time_iter=time_vec;
name = sprintf('%s_time%05s_z%05s','data_1cluster',num2str(time_iter),num2str(z_hight)); 
s  = load(name);
rg = s(1:step:end,2); ripley = s(1:step:end,3); 
ripley=smooth(ripley,2);ripley=ripley*A; 
loglog(rg,ripley,'b<-','linewidth',2,'markersize',10); hold on


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_hight = 7;  
time_iter=time_vec;
name = sprintf('%s_time%05s_z%05s','data_1cluster',num2str(time_iter),num2str(z_hight)); 
s  = load(name);
rg = s(1:step:end,2); ripley = s(1:step:end,3); 
ripley=smooth(ripley,2);ripley=ripley*A; 
loglog(rg,ripley,'ms-','linewidth',2,'markersize',10); hold on


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
z_hight = 10;  
time_iter=time_vec;
name = sprintf('%s_time%05s_z%05s','data_1cluster',num2str(time_iter),num2str(z_hight)); 
s  = load(name);
rg = s(1:step:end,2); ripley = s(1:step:end,3); 
ripley=smooth(ripley,2);ripley=ripley*A; 
loglog(rg,ripley,'gd-','linewidth',2,'markersize',10); hold on

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rg = s(1:step:end,2); 
ripley=pi*rg.*rg;
loglog(rg,ripley,'--','linewidth',4,'markersize',10,'color',[0.5 0.5 0.5]); hold on


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% linear scale
set(gca,'fontsize',24)
    xlabel('$r$', 'Interpreter', 'Latex','Fontsize',24)
    ylabel('$K(r)$', 'Interpreter', 'Latex','Fontsize',24)
grid on
legend('$z=1$','$z=3$','$z=5$','$z=7$','$z=10$','reference line','location','northwest','interpreter', 'latex')
xlim([0 10])
ylim([0 200])


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function s=smooth(s,num)
len = length(s);
iter = num;
for k = 1:iter
s(2:end-1) = (s(1:end-2) + s(3:end))/2;
end


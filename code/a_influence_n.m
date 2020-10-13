%% Peng for Random graph 

%% global setting
T = 200;
% the  number of agents
n = 80;
% select which dataset
data_set = 1;

%% return the regret with convex function and strongly
% convex function  
%  [Reg,RegB,RegTB] = bodyfat_call_n(n, T, 1, data_set);
 
% strongly convex function
[Reg,RegB,RegTB] = bodyfat_call_n(n,T, 2, data_set);

%%  graph part

Reg_ave = zeros(1,n);
RegB_ave = zeros(1,n);
RegTB_ave = zeros(1,n);
for i = 1 : n
    Reg_ave(i) = Reg(i) * 200 / i;
    RegB_ave(i) = RegB(i) * 200 / i;
    RegTB_ave(i) = RegTB(i) * 200 / i;
end
%%
figure 
set(gca, 'FontSize',14);
h1 = line_fewer_markers( 5 : 80, Reg_ave(5:80), 9,'o-b','MarkerSize', 7, 'linewidth', 1.5);
h2 = line_fewer_markers( 5 : 80, RegB_ave(5:80), 9,'d-r','MarkerSize', 7, 'linewidth', 1.5);
h3 = line_fewer_markers( 5 : 80, RegTB_ave(5:80), 9,'s-g','MarkerSize', 7, 'linewidth', 1.5);
% legend([h1 h2 h3],...
%     'Convex: full','Convex: one-point bandit','Convex: two-points bandit');

legend([h1 h2 h3],...
    'Strongly convex: full','Strongly convex: one-point bandit','Strongly convex: two-points bandit');

grid on;
xlabel('N');
zz=ylabel('$$\frac{\max_i~ Reg(i,T)}{N} $$'); 
set(zz,'Interpreter','latex');  


% % strongly convex function
%   [Reg_S1,RegB_S1,RegTB_S1] =  bodyfat_call_probability(n,T, 2, data_set);
%   figure 
% set(gca, 'FontSize',14);
% Reg_S1_ave = zeros(1,9);
% RegB_S1_ave = zeros(1,9);
% RegTB_S1_ave = zeros(1,9);
% for i = 1 : 9
%    Reg_S1_ave(i) = Reg_S1(i) * 200;
%     RegB_S1_ave(i) = RegB_S1(i) * 200;
%     RegTB_S1_ave(i) = RegTB_S1(i)* 200;
% end
% h1 = line_fewer_markers( 0.1:0.1:0.9, Reg_S1_ave, 9,'o-b','MarkerSize', 7, 'linewidth', 1.5);
% h2 = line_fewer_markers( 0.1:0.1:0.9, RegB_S1_ave, 9,'d-r','MarkerSize', 7, 'linewidth', 1.5);
% h3 = line_fewer_markers( 0.1:0.1:0.9, RegTB_S1_ave, 9,'s-g','MarkerSize', 7, 'linewidth', 1.5);
% 
% legend([h1 h2 h3],...
%     'Strongly Convex: full','Strongly Convex: one point bandit','Strongly Convex: two point bandit');
% 
% grid on;
% xlabel('p');
% zz=ylabel('$$\max_i~ Reg(i,T) $$'); 
% set(zz,'Interpreter','latex');  

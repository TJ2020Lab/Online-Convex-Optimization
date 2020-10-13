%% Peng for Random graph 

%% global setting
T = 200;
% the  number of agents
n = 30;
% select which dataset
data_set = 1;

% convex function or strongly convex

%% return the regret with convex function and strongly
% convex function  
[Reg_c1,RegB_c1,RegTB_c1] = bodyfat_call(n, T, 1, data_set);

% strongly convex function
% [Reg_S1,RegB_S1,RegTB_S1] = bodyfat_call(n,T, 2, data_set);

%%  graph part
set(gca, 'FontSize',18);
% h1 = line_fewer_markers( 1:T, Reg_c1, 8,'o-b','MarkerSize', 7, 'linewidth', 1.5);
% h2 = line_fewer_markers( 1:T, RegB_c1, 8,'d-r','MarkerSize', 7, 'linewidth', 1.5);
% h3 = line_fewer_markers( 1:T, RegTB_c1, 8,'s-g','MarkerSize', 7, 'linewidth', 1.5);
h4 = line_fewer_markers( 1:T, Reg_S1, 8,'o--b','MarkerSize', 7, 'linewidth', 1.5);
h5 = line_fewer_markers( 1:T, RegB_S1, 8,'d--r','MarkerSize', 7, 'linewidth', 1.5);
h6 = line_fewer_markers( 1:T, RegTB_S1, 8,'s--g','MarkerSize', 7, 'linewidth', 1.5);

% legend([h1 h2 h3],...
%     'Convex: full','Convex: one point bandit','Convex: two point bandit');
% 
legend([h4 h5 h6],...
    'Strongly convex: full','Strongly convex: one point bandit','Strongly convex: two point bandit');

grid on;
xlabel('Time horizon T');
ylabel('Maximum regret');
% title('dataset = bodyfat; n = 30, p = 0.2')





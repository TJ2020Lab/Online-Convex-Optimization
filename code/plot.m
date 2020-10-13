% %%  graph part
% set(gca, 'FontSize',14);
% h18 = line_fewer_markers( 1:T, Reg_c1_08, 8,'o-b','MarkerSize', 7, 'linewidth', 1.5);
% h15 = line_fewer_markers( 1:T, Reg_c1_05, 8,'o--b','MarkerSize', 7, 'linewidth', 1.5);
% h12 = line_fewer_markers( 1:T, Reg_c1_02, 8,'o-.b','MarkerSize', 7, 'linewidth', 1.5);
% h28 = line_fewer_markers( 1:T, RegB_c1_08, 8,'d-r','MarkerSize', 7, 'linewidth', 1.5);
% h25 = line_fewer_markers( 1:T, RegB_c1_05, 8,'d--r','MarkerSize', 7, 'linewidth', 1.5);
% h22 = line_fewer_markers( 1:T, RegB_c1_02, 8,'d-.r','MarkerSize', 7, 'linewidth', 1.5);
% h38 = line_fewer_markers( 1:T, RegTB_c1_08, 8,'s-g','MarkerSize', 7, 'linewidth', 1.5);
% h35 = line_fewer_markers( 1:T, RegTB_c1_05, 8,'s--g','MarkerSize', 7, 'linewidth', 1.5);
% h32 = line_fewer_markers( 1:T, RegTB_c1_02, 8,'s-.g','MarkerSize', 7, 'linewidth', 1.5);
% % h4 = line_fewer_markers( 1:T, Reg_S1, 8,'o--b','MarkerSize', 7, 'linewidth', 1.5);
% % h5 = line_fewer_markers( 1:T, RegB_S1, 8,'d--r','MarkerSize', 7, 'linewidth', 1.5);
% % h6 = line_fewer_markers( 1:T, RegTB_S1, 8,'s--g','MarkerSize', 7, 'linewidth', 1.5);
% 
% legend([h18 h15 h12 h28 h25 h22 h38 h35 h32],...
%     'Convex: full(p=0.8)','Convex: full(p=0.5)','Convex: full(0.2)',...
%     'Convex: one point bandit(p=0.8)','Convex: one point bandit(p=0.5)','Convex: one point bandit(0.2)',...
%     'Convex: two point bandit(p=0.8)','Convex: two point bandit(p=0.5)','Convex: two point bandit(0.2)');
% 
% grid on;
% xlabel('Time horizon T');
% ylabel('Maximum regret');
%%  graph part
Reg_c1_ave = zeros(1,200);
RegB_c1_ave = zeros(1,200);
RegTB_c1_ave = zeros(1,200);
Reg_S1_ave = zeros(1,200);
RegB_S1_ave = zeros(1,200);
RegTB_S1_ave = zeros(1,200);
% for i = 1 : 9
%     Reg_c1_ave(i) = Reg_c1(i) * 200;
%     RegB_c1_ave(i) = RegB_c1(i) * 200;
%     RegTB_c1_ave(i) = RegTB_c1(i) * 200;
% %     Reg_S1_ave(i) = Reg_S1(i) * i;
% %     RegB_S1_ave(i) = RegB_S1(i) * i;
% %     RegTB_S1_ave(i) = RegTB_S1(i) * i;
% end
% h1 = line_fewer_markers( 1:200, Reg_c1, 8,'o-b','MarkerSize', 7, 'linewidth', 3);
% h2 = line_fewer_markers( 1:200, RegB_c1, 8,'d-r','MarkerSize', 7, 'linewidth', 3);
% h3 = line_fewer_markers( 1:200, RegTB_c1, 8,'s-g','MarkerSize', 7, 'linewidth', 3);
h4 = line_fewer_markers( 1:200, Reg_S1, 8,'o--b','MarkerSize', 7, 'linewidth', 3);
h5 = line_fewer_markers( 1:200, RegB_S1, 8,'d--r','MarkerSize', 7, 'linewidth', 3);
h6 = line_fewer_markers( 1:200, RegTB_S1, 8,'s--g','MarkerSize', 7, 'linewidth', 3);

% legend([h1 h2 h3],...
%     'Convex: full','Convex: one-point bandit','Convex: two-points bandit');
legend([h4 h5 h6],...
    'Strongly convex: full','Strongly convex: one-point bandit','Strongly convex: two-points bandit');

grid on;
xlabel('Time horizon T','FontSize',18);
ylabel('SReg(T)', 'FontSize',18);
set(gca,'FontSize',30);
% title('dataset = bodyfat; n = 30, p = 0.2')



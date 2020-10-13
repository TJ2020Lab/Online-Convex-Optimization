% f1 = imread('n_p05_10.png');
% f2= imread('n_p05_05_s.png');
% f3=imread('d_05.png');
% f4 = imread('d_05_s.png');
% figure(1)
% subplot(2,8,[1,2,9,10])
% % set(gca,'position',[0.5 0.2 1 1])
% imshow(f1);
% subplot(2,8,[3,4,11,12])
% % set(gca,'position',[0.2 0.05 0.5 0.5])
% imshow(f2);
% subplot(2,8,[5,6,13,14])
% % set(gca,'position',[0.4 0.05 0.5 0.5])
% imshow(f3);
% subplot(2,8,[7,8,15,16])
% % set(gca,'position',[0.6 0.05 0.5 0.5])
% imshow(f4);
set(gca, 'FontSize',14);
figure(1)
subplot(1,4,1)

x = 0:75;
       p1 =   0.0001854 ;
       p2 =    -0.03929 ;
       p3 =       2.657;
       p4 =         134;
yy =  p1*x.^3 + p2*x.^2 + p3*x + p4;
h1 = line_fewer_markers( 5 : n, Reg_ave_1(5:n), 9,'o-b','MarkerSize', 7, 'linewidth', 1.5);
h2 = line_fewer_markers( 5 : n, RegB_ave_1(5:n), 9,'d-r','MarkerSize', 7, 'linewidth', 1.5);
h3 = line_fewer_markers( 5 : n, RegTB_ave_1(5:n), 9,'s-g','MarkerSize', 7, 'linewidth', 1.5);
h4 = line_fewer_markers( 5 : n, yy, 9,'d--k','MarkerSize', 7, 'linewidth', 1.5);
legend([h1 h2 h4 h3],...
    'Convex: full','Convex: one-point bandit','Convex: one-point bandit (fit)','Convex: two-points bandit');

% legend([h1 h2 h3],...
%     'Strongly convex: full','Strongly convex: one-point bandit','Strongly convex: two-points bandit');

grid on;
xlabel('N');
zz=ylabel('$$\frac{\max_i~ Reg(i,T)}{N} $$'); 
set(zz,'Interpreter','latex');  



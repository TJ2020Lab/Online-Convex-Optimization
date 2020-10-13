figure 
set(gca, 'FontSize',14);
x = 0:75;
       p1 =   0.0001854 ;
       p2 =    -0.03929 ;
       p3 =       2.657;
       p4 =         134;
yy =  p1*x.^3 + p2*x.^2 + p3*x + p4;
h1 = line_fewer_markers( 5 : n, Reg_ave(5:n), 9,'o-b','MarkerSize', 7, 'linewidth', 1.5);
h2 = line_fewer_markers( 5 : n, RegB_ave(5:n), 9,'d-r','MarkerSize', 7, 'linewidth', 1.5);
h3 = line_fewer_markers( 5 : n, RegTB_ave(5:n), 9,'s-g','MarkerSize', 7, 'linewidth', 1.5);
h4 = line_fewer_markers( 5 : n, yy, 9,'d--k','MarkerSize', 7, 'linewidth', 1.5);
legend([h1 h2 h4 h3],...
    'Convex: full','Convex: one-point bandit','Convex: one-point bandit (fit)','Convex: two-points bandit');

% legend([h1 h2 h3],...
%     'Strongly convex: full','Strongly convex: one-point bandit','Strongly convex: two-points bandit');

grid on;
xlabel('N');
zz=ylabel('$$\frac{\max_i~ Reg(i,T)}{N} $$'); 
set(zz,'Interpreter','latex');  
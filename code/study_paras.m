%%  investigate the paras

function [] = study_paras();
% for non-strongly convex
U = 0.15;
iter = 1000;

N = 1;

c_set = [1/2 3/4];


[a_full,~,~] = long_term_cons(c_set(1),U,iter,1);
[b_full,~,~] = long_term_cons(c_set(2),U,iter,1);

a_avg = zeros(1,iter);
b_avg = zeros(1,iter);
for i = 1 : N
    [a_bandit,~,~] = long_term_cons(c_set(1),U,iter,2);
    [b_bandit,~,~] = long_term_cons(c_set(2),U,iter,2);
    a_avg = a_avg + a_bandit;
    b_avg = b_avg + b_bandit;
end
a_avg = a_avg/N;
b_avg = b_avg/N;

h1 = line_fewer_markers( 1:iter, a_full, 10,'o-b','MarkerSize', 8);
h2 = line_fewer_markers( 1:iter, a_avg, 10,'o--b','MarkerSize', 8);
h3 = line_fewer_markers( 1:iter, b_full, 10,'^-r','MarkerSize', 8);
h4 = line_fewer_markers( 1:iter, b_avg, 10,'^--r','MarkerSize', 8);
legend([h1 h2 h3 h4],'CACV: full(c = 1/2)','CACV: bandit (c = 1/2)',...
    'CACV: full (c = 3/4)','CACV: bandit (c = 3/4)');





%% for strongly convex
U = 0.15;
iter = 1000;
N = 1;


[A_full,~,~] = long_term_cons_SC(U,iter,1,1);
[B_full,~,~] = long_term_cons_SC(U,iter,1,2);

A_avg = zeros(1,iter);
B_avg = zeros(1,iter);

for i = 1 : N
    [A_bandit,~,~] = long_term_cons_SC(U,iter,2,1);
    [B_bandit,~,~] = long_term_cons_SC(U,iter,2,2);
    A_avg = A_avg + A_bandit;
    B_avg = B_avg + B_bandit;
end
A_avg = A_avg/N;
B_avg = B_avg/N;

h1 = line_fewer_markers( 1:iter, A_full, 10,'o-b','MarkerSize', 10);
h2 = line_fewer_markers( 1:iter, A_avg, 10,'o--b','MarkerSize', 10);
h3 = line_fewer_markers( 1:iter, B_full, 10,'^-r','MarkerSize', 10);
h4 = line_fewer_markers( 1:iter, B_avg, 10,'^--r','MarkerSize', 10);
legend([h1 h2 h3 h4],'CACV: full (\rho = 1)','CACV: bandit (\rho = 1)',...
    'CACV: full (\rho = 2)','CACV: bandit (\rho = 2)');






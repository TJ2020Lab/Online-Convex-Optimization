% created by Deming 
%  Revised by Peng for Random graph


function [Reg_T,Reg_TB] = iclr_ltc_REAL_call(c,T)
% read the data
[A,B] = libsvmread('bodyfat_scale.txt');


% parameters
m = 14;
n = 6;
U = .25;
L = -U;
T = 200;
c = 1/2;
b = c/3;
rho = 1;

c_eta_reg = 20;  % regularization
c_beta_step = 30;  % full step
c_beta_step_B = 20;  % bandit step
c_epsi = 1;  % exploration
eta_reg = zeros(1,T);
beta_step = zeros(1,T);
beta_step_B = zeros(1,T);
epsi = zeros(1,T);
% pi = epsi/R;

for s = 1 : T
    eta_reg(s) = 1/(c_eta_reg*(s^c));
    beta_step(s) = 1/(c_beta_step*(s^c));
    beta_step_B(s) = 1/(c_beta_step_B*(s^c));
    epsi(s) = 1/(c_epsi*(s^b));
end

f_opt = zeros(T,1);
pro1 = zeros(m,m);
pro2 = zeros(m,1);
pro3 = 0;
% store all the samples
B_randperm = zeros(T,m,n);
A_randperm = zeros(T,n);
for t = 1 : T
    for p = 1 : n
        t_fake = randi([1 252],1);      
        pro1 = pro1 + B(t_fake,:)'*B(t_fake,:);
        pro2 = pro2 + A(t_fake)*B(t_fake,:)';
        pro3 = pro3 + A(t_fake)^2;    
        B_randperm(t,:,p) = B(t_fake,:);
        A_randperm(t,p) = A(t_fake);
    end
    x_opt = (pro1 + 2*rho*n*t*eye(m))\pro2;
    x_opt = box_proj(m,x_opt,L,U);
    f_opt(t) = (1/2)*x_opt'*pro1*x_opt - pro2'*x_opt + (1/2)*pro3 +...
        n*t*rho*norm(x_opt,2)^2;
end



% compute the regrets
Reg_T = zeros(1,T);
Reg_TB = zeros(1,T);
% cacvs
CACV_T = zeros(1,T);
CACV_TB = zeros(1,T);
% docg regrets
Reg_T_docg = zeros(1,T);
Reg_TB_docg = zeros(1,T);

for t = 1 : T
    [Reg_T(t),CACV_T(t),Reg_T_docg(t)] =  iclr_ltc_REAL(U,t,rho,...
        B_randperm(1:t,:,:),A_randperm(1:t,:),f_opt(t),...
        beta_step(t),beta_step_B(t),epsi(t),eta_reg(t),1);
    [Reg_TB(t),CACV_TB(t),Reg_TB_docg(t)] = iclr_ltc_REAL(U,t,rho,...
        B_randperm(1:t,:,:),A_randperm(1:t,:),f_opt(t),...
        beta_step(t),beta_step_B(t),epsi(t),eta_reg(t),2);
end


set(gca, 'FontSize',14);
h1 = line_fewer_markers( 1:T, Reg_T, 10,'o-b','MarkerSize', 10, 'linewidth', 2);
h2 = line_fewer_markers( 1:T, Reg_TB, 10,'d-r','MarkerSize', 10, 'linewidth', 2);
h3 = line_fewer_markers( 1:T, Reg_T_docg, 10,'o-m','MarkerSize', 10, 'linewidth', 2);
h4 = line_fewer_markers( 1:T, Reg_TB_docg, 10,'d-c','MarkerSize', 10, 'linewidth', 2);
h5 = line_fewer_markers( 1:T, CACV_T, 10,'o--b','MarkerSize', 10, 'linewidth', 2);
h6 = line_fewer_markers( 1:T, CACV_TB, 10,'d--r','MarkerSize', 10, 'linewidth', 2);

legend([h1 h2 h3 h4 h5 h6],...
    'DOCO: maximum regret (full)','DOCO: maximum regret (bandit)',...
    'D-OCG: maximum regret (full)','D-OCG: maximum regret (bandit)',...
    'DOCO: CACV (full)','DOCO: CACV (bandit)');

grid on;
xlabel('Time horizon T');
% ylabel('Maximum regret');
title('dataset = bodyfat; \rho = 1')
set(gca, 'YScale', 'linear');







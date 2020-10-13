%% 

function [Reg_T,Reg_TB,CACV_T,CACV_TB] = iclr_ltc_call(c,T)
n = 6;
m = 4;
U = 0.5;
b = c/3;

c_eta_reg = 1;  % regularization
c_beta_step = .5;  % full and bandit use the same step
c_epsi = 1;  % exploration
eta_reg = zeros(1,T);
beta_step = zeros(1,T);
epsi = zeros(1,T);

for s = 1 : T
    eta_reg(s) = 1/(c_eta_reg*(s^c));
    beta_step(s) = 1/(c_beta_step*(s^c));
    epsi(s) = 1/(c_epsi*(s^b));
end

% the input data
rng(1,'twister');
pro_a = 1 + (-2)*rand(m,n,T);

% the output
x_zero = zeros(m,1);
for i = 1 : m/2
    x_zero(i) = 1;
end

rng(1,'twister');
var_noise = randn(n,T); % time-varying noise
b_out = zeros(n,T);
for j = 1 : n
    for r = 1 : T
        b_out(j,r) = pro_a(:,j,r)'*(x_zero) + var_noise(j,r);
    end
end

f_opt = zeros(T,1);
pro1 = zeros(m,m);
pro2 = zeros(m,1);
pro3 = 0;
for t = 1 : T
    for p = 1 : n
        pro1 = pro1 + pro_a(:,p,t)*pro_a(:,p,t)';
        pro2 = pro2 + b_out(p,t)*pro_a(:,p,t);
        pro3 = pro3 + b_out(p,t)^2;
    end
    x_opt = pro1\pro2;
    f_opt(t) = (1/2)*x_opt'*pro1*x_opt - pro2'*x_opt + (1/2)*pro3;
end

% compute the regrets
Reg_T = zeros(1,T);
Reg_TB = zeros(1,T);
% cacvs
CACV_T = zeros(1,T);
CACV_TB = zeros(1,T);

for t = 1 : T
    [Reg_T(t),CACV_T(t)] =  iclr_ltc(m,n,U,t,...
        pro_a(:,:,1:t),b_out(:,1:t),f_opt(t),...
        beta_step(t),epsi(t),eta_reg(t),1);
    [Reg_TB(t),CACV_TB(t)] = iclr_ltc(m,n,U,t,...
        pro_a(:,:,1:t),b_out(:,1:t),f_opt(t),...
        beta_step(t),epsi(t),eta_reg(t),2);
end






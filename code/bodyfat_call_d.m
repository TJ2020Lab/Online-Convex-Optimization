% created by Deming 
%  Revised by Peng for Random graph


function [Reg_F_ave,Reg_B_ave,Reg_TB_ave] = bodyfat_call_d(n,T,flag_c)
d = 5:5:100;
Reg_F_ave = zeros(1, length(d));
Reg_B_ave = zeros(1, length(d));
Reg_TB_ave = zeros(1, length(d));
if flag_c==1
    % convex function
    rho=0;
    % the step-size parameter
    c=1/2;
else
    % strongly convex function
    rho = 1;
    % the step-size parameter
    c=1;
end

for d_i = 1:length(d)
    m = d(d_i);
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


% parameters
R = 5;
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
    x_opt = (pro1 + 2*rho*n*t*eye(m))\pro2;
    f_opt(t) = (1/2)*x_opt'*pro1*x_opt - pro2'*x_opt + (1/2)*pro3+...
        n*t*rho*norm(x_opt,2)^2;
end



delta = zeros(2, T);
% bandit case parameter
for t = 1 : T
    if flag_c==1
        % convex function
        delta(1,t) = 1/(t^0.25);
        delta(2,t) = 1/(t^0.5);
    else
        % strongly convex function
        delta(1,t) = ((1+log(t))/t)^(1/3);
        delta(2,t) = log(t)/t;
    end
end

xi = delta / R;

% pi = epsi/R;


%% Regret vs d

prob = 0.5;

% k-regular graph with k=3
kr=3;
ad_k = kregular(n, kr);

Reg_F_sum = 0;
Reg_B_sum = 0;
Reg_TB_sum = 0;
for ave_i = 1 : 20
        % ER-Random graph
        ad = zeros(n, n, T);
        wm = zeros(n, n, T);
        
%         % complete
%         for i = 1 : T
%             ad(:,:,i) = randomGraph(n, prob);
%         end
%         
%         % star
%         for i = 1 : T
%             for j = 2 : n
%                 if rand < prob
%                     ad(1,j,i) = 1;
%                     ad(j,1,i) = ad(1,j,i);
%                 end
%             end
%         end
%        
%         for t = 1 : T
%             for i = 1 : n
%                 for j = i+1 : n
%                     if ad(i,j,t) == 1
%                         wm(i,j,t) = 1 / n;
%                         wm(j,i,t) = wm(i,j,t);
%                     end
%                 end
%                 wm(i,i,t) = 1 - sum(wm(i,:,t));
%             end
%         end
         

        % k-regular, k=2 is the cycle graph
        for i = 1 : T
            for j = 1 : n
                for l = j+1 : n
                    if ad_k(j,l) == 1
                        if rand < prob
                            ad(j,l,i) = 1;
                            ad(l,j,i) = ad(j,l,i);
                        end
                    end
                end
            end
        end
      %  wm of k-regular
        for t = 1 : T
            for i = 1 : n
                for j = i+1 : n
                    if ad(i,j,t) == 1
                        wm(i,j,t) = 1 / (kr+1);
                        wm(j,i,t) = wm(i,j,t);
                    end
                end
                wm(i,i,t) = 1 - sum(wm(i,:,t));
            end
        end



        % compute the regrets
        t = T;
        [Reg_F] =  iclr_ltc_REAL_d(m,n,t,rho,...
            pro_a(:,:,1:t),b_out(:,1:t),f_opt(t),...
            c,flag_c, 0,0,0,wm(:,:,1:t),R);
        [Reg_B] = iclr_ltc_REAL_d(m,n,t,rho,...
            pro_a(:,:,1:t),b_out(:,1:t),f_opt(t),...
            c,flag_c,delta(1,t),xi(1,t),1, wm(:,:,1:t),R);
        [Reg_TB] = iclr_ltc_REAL_d(m,n,t,rho,...
            pro_a(:,:,1:t),b_out(:,1:t),f_opt(t),...
            c,flag_c,delta(2,t),xi(1,t),2, wm(:,:,1:t),R);
        Reg_F_sum = Reg_F_sum + Reg_F;
         Reg_B_sum = Reg_B_sum + Reg_B;
         Reg_TB_sum = Reg_TB_sum + Reg_TB;
end
        Reg_F_ave(d_i) =  Reg_F_sum / 20;
        Reg_B_ave(d_i) = Reg_B_sum / 20;
        Reg_TB_ave(d_i) =  Reg_TB_sum / 20;
end
end

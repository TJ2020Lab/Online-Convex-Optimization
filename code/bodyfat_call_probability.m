% created by Deming 
%  Revised by Peng for Random graph


function [Reg_F_ave,Reg_B_ave,Reg_TB_ave] = bodyfat_call_probability(n,T,flag_c,data_set)
% read the data
if data_set==1
[A,B] = libsvmread('bodyfat_scale.txt');
else 
[A,B] = libsvmread('abalone_scale.txt'); 
end

[Sample_number, feature_dim]=size(B);

% parameters
R = 5;
% The number of features, the dimension of x 
m = feature_dim;
samples = Sample_number;

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

f_opt = zeros(T,1);
pro1 = zeros(m,m);
pro2 = zeros(m,1);
pro3 = 0;
% store all the samples
B_randperm = zeros(T,m,n);
A_randperm = zeros(T,n);
for t = 1 : T
    for p = 1 : n
        t_fake = randi([1 samples],1);      
        pro1 = pro1 + B(t_fake,:)'*B(t_fake,:);
        pro2 = pro2 + A(t_fake)*B(t_fake,:)';
        pro3 = pro3 + A(t_fake)^2;    
        B_randperm(t,:,p) = B(t_fake,:);
        A_randperm(t,p) = A(t_fake);
    end
    x_opt = (pro1 + 2*rho*n*t*eye(m))\pro2;
    f_opt(t) = (1/2)*x_opt'*pro1*x_opt - pro2'*x_opt + (1/2)*pro3 +...
        n*t*rho*norm(x_opt,2)^2;
end
%% Regret vs probability
Reg_F_ave = zeros(1, 9);
Reg_B_ave = zeros(1, 9);
Reg_TB_ave = zeros(1, 9);


% k-regular graph with k=5
kr=5;
ad_k = kregular(n, kr);


for prob = 0.1: 0.1 : 0.9
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
        [Reg_F] =  iclr_ltc_REAL(m,n,t,rho,...
            B_randperm(1:t,:,:),A_randperm(1:t,:),f_opt(t),...
            c,flag_c, 0,0,0,wm(:,:,1:t),R);
        [Reg_B] = iclr_ltc_REAL(m,n,t,rho,...
            B_randperm(1:t,:,:),A_randperm(1:t,:),f_opt(t),...
            c,flag_c,delta(1,t),xi(1,t),1, wm(:,:,1:t),R);
        [Reg_TB] = iclr_ltc_REAL(m,n,t,rho,...
            B_randperm(1:t,:,:),A_randperm(1:t,:),f_opt(t),...
            c,flag_c,delta(2,t),xi(1,t),2, wm(:,:,1:t),R);
        ind = fix(prob * 10);
        Reg_F_ave(ind) = Reg_F_ave(ind) + Reg_F;
        Reg_B_ave(ind) = Reg_B_ave(ind) + Reg_B;
        Reg_TB_ave(ind) = Reg_TB_ave(ind) + Reg_TB;
    end
end
Reg_F_ave = Reg_F_ave / 20;
Reg_B_ave = Reg_B_ave / 20;
Reg_TB_ave = Reg_TB_ave / 20;
end








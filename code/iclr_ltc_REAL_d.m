%% 

function [reg_max] = iclr_ltc_REAL_d(m,n,iter,rho,...
    D_feature,D_label,f_obj, c,flag_c,delta, xi, full_or_bandit, WM,R)
% parameters

% initial decisions
rng(1,'twister');
x_dec = rand(m,n);
reg_vec = zeros(n,1);

f_val = zeros(n,1);
for i = 1 : n
    for j = 1 : n
        f_val(i) = f_val(i) +...
            (1/2)*( D_feature(:,j,1)'*x_dec(:,i) - D_label(j,1) )^2+...
                        rho*norm(x_dec(:,i),2)^2;
    end
end



% convex
c_eta_step = 0.5;  % full
c_eta_B_step = 0.5;  % one point bandit
c_eta_TB_step = 0.5;  % two point bandit
% strongly convex
s_eta_step = 0.5;  % full
s_eta_B_step = 0.5;  % one point bandit
s_eta_TB_step = 0.5;  % two point bandit

eta_step = zeros(1,iter);
eta_B_step = zeros(1,iter);
eta_TB_step = zeros(1,iter);

for s = 1 : iter
    if flag_c==1
        eta_step(s) = 1/(c_eta_step*(s^c));
        eta_B_step(s) = delta /(c_eta_B_step*(s^c));
        eta_TB_step(s) = 1/(c_eta_TB_step*(s^c));
    else
        eta_step(s) = 1/(s_eta_step*(s^c));
        eta_B_step(s) = 1/(s_eta_B_step*(s^c));
        eta_TB_step(s) = 1/(s_eta_TB_step*(s^c));
    end 
end



for t = 2 : iter
    WM_t = WM(:,:,t-1);
    for j = 1 : n
        % compute the subgradient
        subg = ( D_feature(:,j,t)'*x_dec(:,j) - D_label(j,t))*D_feature(:,j,t)+...
            2*rho*x_dec(:,j);
        
        % compute the ONE-point feedback
        ran_vec = randn(m,1);
        ran_vec = ran_vec./norm(ran_vec,2);
        x_dec_t1 = x_dec(:,j) + delta*ran_vec;
        subg_ban1 = ( (1/2)*( D_feature(:,j,t)'*x_dec_t1 - D_label(j,t) )^2 +...
            rho*norm(x_dec_t1,2)^2 )*(m/(delta))*ran_vec;
        
        % compute the TWO-point feedback
        x_dec_t2 = x_dec(:,j) - delta*ran_vec;
        subg_ban2 = (( (1/2)*( D_feature(:,j,t)'*x_dec_t1 - D_label(j,t) )^2 +...
            rho*norm(x_dec_t1,2)^2 ) - ((1/2)*( D_feature(:,j,t)'*x_dec_t2 - D_label(j,t) )^2 +...
            rho*norm(x_dec_t2,2)^2 ))*(m/(2*delta))*ran_vec;        
        

        
        % full_or_bandit == 0, run full;
        switch(full_or_bandit)
            case 0
                % DOCO - full
                x_dec(:,j) = x_dec(:,j) - eta_step(t) * subg;
            case 1
                % DOCO - one bandit
                x_dec(:,j) = x_dec(:,j) - eta_B_step(t) * subg_ban1;
            case 2
                % DOCO - two bandit
                x_dec(:,j) = x_dec(:,j) - eta_TB_step(t) * subg_ban2;
            %             R = (1-pi)*R;
        end
    end
    
    % average consensus
    x_dec = x_dec*WM_t';
    
    % projection
    for l = 1 : n
        x_dec(:,l) = ( (1-xi)*R / max(norm(x_dec(:,l)),(1-xi)*R) )*x_dec(:,l);
    end
    
    
    % sum all the function values
    for i = 1 : n
        for p = 1 : n
            f_val(i) = f_val(i) +...
                (1/2)*( D_feature(:,p,t)'*x_dec(:,i) - D_label(p,t) )^2+...
                rho*norm(x_dec(:,i),2)^2;
           
        end
    end
    
end


% calculate the regrets for every node i
for i = 1 : n
    reg_vec(i) = (f_val(i) - f_obj);
end
reg_max = max(reg_vec);
reg_max = reg_max / iter;


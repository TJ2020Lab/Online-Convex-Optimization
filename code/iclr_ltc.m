%%
function [reg_max,cacv] = iclr_ltc(m,n,U,iter,...
    pro_a,b_out,f_obj,...
    beta_step,epsi,eta_reg,...
    full_or_bandit)

% parameters
p = 2*m;           % number of inequalities
L = -U;
R = sqrt(m)*U;

%********************************************************************
% a time-varying network
%********************************************************************
WM = zeros(6,6,4);
WM(:,:,1) = [1/2 0 0 1/2 0 0;
    0 1/2 0 0 1/2 0;
    0 0 1/2 0 0 1/2;
    1/2 0 0 1/2 0 0;
    0 1/2 0 0 1/2 0;
    0 0 1/2 0 0 1/2];

WM(:,:,2) = [2/3 1/3 0 0 0 0;
    1/3 1/3 1/3 0 0 0;
    0 1/3 2/3 0 0 0;
    0 0 0 2/3 1/3 0;
    0 0 0 1/3 1/3 1/3;
    0 0 0 0 1/3 2/3];

WM(:,:,3) = [2/3 1/3 0 0 0 0;
    1/3 1/3 1/3 0 0 0;
    0 1/3 1/3 0 0 1/3;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 1/3 0 0 2/3];

WM(:,:,4) = [2/3 0 0 1/3 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    1/3 0 0 1/3 1/3 0;
    0 0 0 1/3 1/3 1/3;
    0 0 0 0 1/3 2/3];


% initial decisions
rng(1,'twister');
x_dec = rand(m,n);
reg_vec = zeros(n,1);            % regret vectors for all the nodes


dual_dec = zeros(p,n);
% calculate the cacv for the 1st round
for y = 1 : n
    for s = 1 : m
        cacv = max(L - x_dec(s,y),0) + max(x_dec(s,y) - U,0);
    end
end


f_val = zeros(n,1);
for i = 1 : n
    for j = 1 : n
        f_val(i) = f_val(i) +...
            (1/2)*( pro_a(:,j,1)'*x_dec(:,i) - b_out(j,1) )^2;
    end
end


for t = 2 : iter
    for j = 1 : n
        % compute the subgradient
        subg = ( pro_a(:,j,t)'*x_dec(:,j) - b_out(j,t))*pro_a(:,j,t);
        
        % compute the ONE-point feedback
        ran_vec = randn(m,1);
        ran_vec = ran_vec./norm(ran_vec,2);
        x_dec_t1 = x_dec(:,j) + epsi*ran_vec;
        subg_ban = (1/2)*( pro_a(:,j,t)'*x_dec_t1 - b_out(j,t) )^2*...
            (m/(epsi))*ran_vec;
        
        % full_or_bandit == 1, run full;
        if full_or_bandit == 1
            x_dec(:,j) = x_dec(:,j) - beta_step*( subg +...
                dual_dec(:,j)'*iclr_sub_ineqs(m,x_dec(:,j),L,U) );
        else
            x_dec(:,j) = x_dec(:,j) - beta_step*( subg_ban +...
                dual_dec(:,j)'*iclr_sub_ineqs(m,x_dec(:,j),L,U) );     
            %             R = (1-pi)*R;
        end    
    end
    
    % average consensus
    x_dec = x_dec*WM(:,:,rem(t-1,4)+1)';
    
    % projection and dual update
    for l = 1 : n
        x_dec(:,l) = ( R/max(norm(x_dec(:,l)),R) )*x_dec(:,l);
        for s = 1 : m
            dual_dec(s,l) = max(L - x_dec(s,l),0)/eta_reg;
            dual_dec(s+m,l) = max(x_dec(s,l)-U,0)/eta_reg;
        end
    end
    
    % calculate the cacv for every time t
    for y = 1 : n
        for u = 1 : m
            cacv = cacv + max(L - x_dec(u,y),0) + max(x_dec(u,y) - U,0);
        end
    end
    
    % sum all the function values
    for i = 1 : n
        for p = 1 : n
            f_val(i) = f_val(i) +...
                (1/2)*( pro_a(:,p,t)'*x_dec(:,i) - b_out(p,t) )^2;
        end
    end
    
end

% calculate the regrets for every node i
for i = 1 : n
    reg_vec(i) = f_val(i) - f_obj;
end
reg_max = max(reg_vec);





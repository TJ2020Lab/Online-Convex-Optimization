
%% created by Deming 

function [cacv_vec,reg_max,reg_min] = ...
    iclr_ltc_SC(U,iter,full_or_bandit,rho);


%********************************************************************
% problem data
%********************************************************************
n = 6;             % number of nodes
m = 4;            % dimension
% iter = 50;         % number of iterations
p = 2*m;           % number of inequalities
% U = .5;
L = -U;
R = sqrt(m)*U;
% rho = 1;           % strongly convex para.



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

%********************************************************************
% the time-varying network
%********************************************************************
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


%********************************************************************
% the algorithm
%********************************************************************
% the input data
rng(1,'twister');
pro_a = 1 + (-2)*rand(m,n,iter);
% pro_a = randn(m,n,iter);   % normal distribution

% the output data
x_zero = zeros(m,1);
for i = 1 : m/2
    x_zero(i) = 1;
end
rng(1,'twister');
var_noise = randn(n,iter); % time-varying noise
b_out = zeros(n,iter);
for j = 1 : n
    for r = 1 : iter
        b_out(j,r) = pro_a(:,j,r)'*x_zero + var_noise(j,r);
    end
end



% initial decisions
x_dec = zeros(m,n);
dual_dec = zeros(p,n);
x_dec_all = zeros(m,n,iter);         % all the decisions
reg_vec = zeros(n,iter);             % regret vectors
cacv_vec = zeros(1,iter);            % cacv vectors
f_temp_all = zeros(n,1);

for p = 1 : n
    for h = 1 : n
        f_temp_all(p) = f_temp_all(p)...
            + (1/2)*( pro_a(:,h,1)'*x_dec_all(:,p,1) - b_out(h,1) )^2 +...
            rho*norm(x_dec_all(:,p,1),2)^2;
    end
end



% parameters in Bandit setting
b = 1/3;
ep = 2;
epsi = 1/(ep*iter^b);
% pi = epsi/R;


for t = 2 : iter
    %     tic;
    beta_step = 1/(rho*t);
    eta_reg = 1/(rho*t);
    
    p_dec = zeros(m,n);
    for j = 1 : n
        % compute the subgradient
        subg = ( pro_a(:,j,t)'*x_dec(:,j) - b_out(j,t))*pro_a(:,j,t) +...
            2*rho*x_dec(:,j);
        
        % compute the ONE-point feedback
        ran_vec = randn(m,1);
        ran_vec = ran_vec./norm(ran_vec,2);
        x_dec_t1 = x_dec(:,j) + epsi*ran_vec;
        subg_ban = ( (1/2)*( pro_a(:,j,t)'*x_dec_t1 - b_out(j,t) )^2 +...
           rho*norm(x_dec_t1,2)^2 )*(m/(epsi))*ran_vec;
        
        % full_or_bandit == 1, run full;
        if full_or_bandit == 1
            
            p_dec(:,j) = x_dec(:,j) - beta_step*( subg +...
                dual_dec(:,j)'*sub_ineqs(m,x_dec(:,j),L,U) );
            
            % full_or_bandit \= 1, run bandit;
        else
            
            p_dec(:,j) = x_dec(:,j) - beta_step*( subg_ban +...
                dual_dec(:,j)'*sub_ineqs(m,x_dec(:,j),L,U) );
            
%             R = (1-pi)*R;
            
        end
        
    end
    
    % average consensus
    p_dec = p_dec*WM(:,:,rem(t-1,4)+1)';

    
    for l = 1 : n
        x_dec(:,l) = ( R/max(norm(p_dec(:,l)),R) )*p_dec(:,l);
        for s = 1 : m
            dual_dec(s,l) = max(L - x_dec(s,l),0)/eta_reg;
            dual_dec(s+m,l) = max(x_dec(s,l)-U,0)/eta_reg;
        end
    end
    
    % calculate the cacv for every time t
    for y = 1 : n
        for u = 1 : m
            cacv_vec(t) = cacv_vec(t-1) + max(L - x_dec(u,y),0) +...
                max(x_dec(u,y) - U,0);
        end
    end
    cacv_vec(t) = cacv_vec(t);
    
    
    x_dec_all(:,:,t) = x_dec;
    
    

    pro_1 = zeros(m,m);
    pro_2 = zeros(m,1);
    pro_3 = 0;
    
    for d = 1 : t
        for p = 1 : n
            pro_1 = pro_1 + pro_a(:,p,d)*pro_a(:,p,d)';
            pro_2 = pro_2 + b_out(p,d)*pro_a(:,p,d);
            pro_3 = pro_3 + b_out(p,d)^2;
        end
    end
    

    x_min = inv( 2*rho*n*t*eye(m) +  pro_1 )*pro_2;
    x_min = box_proj(m,x_min,L,U);
    f_obj = (1/2)*x_min'*pro_1*x_min - pro_2'*x_min + (1/2)*pro_3 +...
        n*t*rho*norm(x_min,2)^2;
    
 
    
    %     tic;
    
    
    % calculate the regrets for every node p
    for p = 1 : n
        for h = 1 : n
            f_temp_all(p) = f_temp_all(p)...
                + (1/2)*( pro_a(:,h,t)'*x_dec_all(:,p,t) - b_out(h,t) )^2 +...
                rho*norm(x_dec_all(:,p,t),2)^2;
        end
        
        reg_vec(p,t) = f_temp_all(p) - f_obj ;   % the time average regret
        
    end
    
    
    
    %     elapsed_3 = toc;
    
end


%********************************************************************
% compute the first regret
%********************************************************************
pro_1_1st = zeros(m,m);
pro_2_1st = zeros(m,1);
pro_3_1st = 0;

for p = 1 : n
    pro_1_1st = pro_1_1st + pro_a(:,p,1)*pro_a(:,p,1)';
    pro_2_1st = pro_2_1st + b_out(p,1)*pro_a(:,p,1);
    pro_3_1st = pro_3_1st + b_out(p,1)^2;
end

cvx_begin
variable x_min_1(m);
expression f_obj_1(1);

f_obj_1 = (1/2)*x_min_1'*pro_1_1st*x_min_1 - pro_2_1st'*x_min_1 +...
    (1/2)*pro_3_1st + n*t*rho*x_min_1'*x_min_1;

minimize( f_obj_1 );
cvx_end

x_min_1 = box_proj(m,x_min_1,L,U);
f_obj_1 = (1/2)*x_min_1'*pro_1_1st*x_min_1 - pro_2_1st'*x_min_1 + (1/2)*pro_3_1st +...
    n*rho*norm(x_min_1,2)^2;


for p = 1 : n
    f_1 = 0;
    for h = 1 : n
        f_1 = f_1 + (1/2)*( pro_a(:,h,1)'*x_dec_all(:,p,1) - b_out(h,1) )^2 +...
            rho*norm(x_dec_all(:,p,1),2)^2;
    end
    reg_vec(p,1) = f_1 - f_obj_1;
end


reg_max =  max(reg_vec);
reg_min =  min(reg_vec);







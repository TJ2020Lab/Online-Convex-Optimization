%% test for distributed resource allocation algorithm
%% test for stochastic resource allocation

%   min \sum_{i=1}^N  f_i(x_i)=E_{\omega_i}x_i(Q_i+\Omega_i)x_i +
%   (c_i+\theta_i)x_i
%   \sum x_i = \sum d_i+\delta_i

clear all
clc
clf 

%% initial parameter setting
%% number of the agents
N=10;

% The probability for the random graph 
pmin=0.05;
pmax=0.1;
% dimension of the decision vector
m=3;

%The algorithm iteration times lens
len=2000;

tic
% test number

tt=1000;
% simulation results
% optimal norm     ||P(S)|| 
op_C=zeros(1, tt);
% resource balance || 1^T_n\otimes I_m (X-D) ||
b_C=zeros(1, tt);  
%censensus result  ||L \Lambda ||
cen_C=zeros(1,tt);
% the eigenvalue of the average grpah s_2(\bar{L})
se_R=zeros(1,tt); 

opsample  =  zeros(len,tt);  % optimal norm      ||  P(S) ||
bsample   =  zeros(len,tt);  %  resource balance || 1_n^T\otimes I_m (X-D)  ||
censample =  zeros(len,tt);  % The censensus erro for lambda || \bar{L} \Lambda  ||


% for flag whether this round is success
tesuc=zeros(tt,1);

for tnm=1:1:tt

%% Generate the cost functions 
% f_i(x_i)= x_i Q_i x_i +  c_i x_i 
% Q = [Q_1,...,Q_N]  \in   R^{m}*{mN}
% C =  [c_1;c_2;c_3;...;c_N] \in R^{mN}*1
% D =  [d_1;d_2;....;d_N] \in R^{mN}*1

Q = zeros(m,m*N);
C = zeros(m*N,1);
D = zeros(m*N,1);

% generation the data Q_i,c_i for each agent i
% Generate the resource vector
for i=1:N
     % for Q
     randd=2+7.*rand;
     tempd = randd*diag(rand(m,1));
     Ut = orth(rand(m,m));
     temp1 = Ut' * tempd * Ut;
     Q(:,(i-1)*m+1:i*m)=temp1;
     % for C
     randd=4+5.*rand;
     temp2=randd*rand(m,1);
     C((i-1)*m+1:i*m,1)=temp2;
     % for D
     randd=3+2*rand;
     temp3=randd.*rand(m,1);
     D((i-1)*m+1:i*m,1)=temp3;
         
end
   
%% generate NUMG random graph 
% random chose one during the algorithm 
% make sure that the average graph is check

NUMG=30;

% store all the Laplacian matrix 
LT   = zeros(N,NUMG*N);
ADJT = zeros(N,NUMG*N);

% the average Laplacian
ADJS=zeros(N,N);

Lavg=zeros(N,N);

  flag=0; % For correctly generat 10 graph 
  while flag==0
      for i=1:NUMG
       % the graph edge probability is random
       p  = rand*(pmax-pmin)+pmin;
      adj = randomGraph(N,p);
      
      Ltemp = laplacianMatrix(adj);
      
      LT(:,(i-1)*N+1:i*N)=Ltemp;
      ADJT(:,(i-1)*N+1:i*N)=adj;
      
      Lavg = Lavg + Ltemp;
      ADJS = ADJS + adj;
      end
      ADJS=1 / NUMG * ADJS;
      Lavg=1 / NUMG * Lavg;
      LADJS=laplacianMatrix(ADJS);
      if (norm(LADJS-Lavg)<0.0001)
         disp('Correctly generation the average laplacian')
      end
      if isConnected(ADJS)==1
          flag=1;
      else 
          flag=0;
      end            
  end
  
 % whether the average graph is connected
 egn= graphSpectrum(ADJS);
 %disp eignvalue
 disp('second eigen')
 egn2=egn(N-1,1)
 if egn(N-1,1)> 0+eps
     disp('Average graph connected ! Right')
 else
     disp('Wrong communication graph!')
 end
 se_R(1,tnm)=egn2;
 
%% algorithm part
 %initialization
 x   =  zeros(m*N, len+1);
 lam =  zeros(m*N, len+1);
 z   =  zeros(m*N, len+1);
 
 % observation model for delta_i
 %
 x(:,1)=D + 0.5*random('norm',0,1,m*N,1);
 % algorithm part

 for k=1:1:len
     % step-size
     
     alphak=2/(k+1)^0.95;
     
     %gradient
     % generation the gradient perturbation noise 
     % Generation the noise observation for local resource
     QN=random('norm',0,0.5, m, m*N);
     CN=random('norm',0,0.5, m*N,1);
     DN=random('norm',0,0.5, m*N,1);
     
     % generate the noise observation for the gradient information
     QP=Q + QN;
     CP=C + CN;
     DP=D + DN;
     
     
     %iterations for x
     gg=zeros(m*N,1);
     for ti=1:N
     gg( (ti-1)*m+1:ti*m,1) = QP(:,(ti-1)*m+1:ti*m)*x( (ti-1)*m+1:ti*m,k);
     end
     gg=gg + CP;
     
     x(:,k+1)   =  x(:,k)   + alphak.* ( - gg + lam(:,k) ); 
     
     % select one random graph from LT
     RGS   =  ceil(rand*NUMG);
     LR    =  LT(:,(RGS-1)*N+1:RGS*N);
     adjR  =  ADJT(:,(RGS-1)*N+1:RGS*N);
     
     % generate the random communication noise 
     
     zeta=zeros(m*N,1);
     
     epln=zeros(m*N,1);
     
     for ti=1:N
         noise1=random('norm',0,1, m*N,1);
         noise2=random('norm',0,1, m*N,1);
         
         zeta((ti-1)*m+1:ti*m,1) = kron(adjR(ti,:),eye(m))*noise1;
         epln((ti-1)*m+1:ti*m,1) = kron(adjR(ti,:),eye(m))*noise2;
     end
 
     lam(:,k+1)  =  lam(:,k)  +  alphak.*( - kron(LR,eye(m)) * (lam(:,k)+ z(:,k))- zeta - epln + (DP-x(:,k)));
     z(:,k+1)    =  z(:,k)    +  alphak.*( kron(LR,eye(m))* lam(:,k)+zeta);
     
      
     % store some performance criterion
     
     % || L\Lambda ||
     cenerr=kron(Lavg,eye(m))*lam(:,k);
     censample(k,tnm) = norm(cenerr);
     
     % || P(S)|| 
     gg=zeros(m*N,1);
     for ti=1:N
     gg( (ti-1)*m+1:ti*m,1) = Q(:,(ti-1)*m+1:ti*m)*x( (ti-1)*m+1:ti*m,k);
     end
     gg=gg + C;
     
     tempop1 =  -gg+lam(:,k);
     tempop2 =  -kron(Lavg,eye(m))*(z(:,k))-cenerr + (D-x(:,k));
     tempop3 =  cenerr;
     
     opsample(k,tnm)=norm([tempop1;tempop2;tempop3]);
     
     
     %  || x - d ||
     
     bsample(k,tnm)= norm(kron(ones(1,N),eye(m))*(x(:,k)-D));
     
     
     
  end
 
 
 %% test the result for this round 
 if (bsample(len,tnm)+opsample(len,tnm)+censample(len,tnm))<20
     tesuc(tnm,1)=1;
 else
     tesuc(tnm,1)=0;
 end
 
 % store the final result for this round
 % optimal norm     ||P(S)|| 
op_C(1,tnm)=opsample(len, tnm);
% resource balance || X-D ||
b_C(1,tnm)=bsample(len, tnm);  
%censensus result  ||L \Lambda ||
cen_C(1,tnm)=censample(len,tnm);
 
end
 

toc 



%% The graph part show the final result

for k=1:len
    t(k,1)=k;
end


% figure(2)
% plot(t(1:len,1),consensusnorm(1:len,1),'g-','LineWidth',4.5)
% hold on
% xlabel('Time');
% title('The trajectories of || L\Lambda||_2');
% 
% legend('||L\Lambda ||_2')
% % %ylabel('The consensus error ||L\Lambda||_2^2 of five agents')
% 
% ll=100;
% tl=t(len-ll:len,1);
% cl=consensusnorm(len-ll:len,1);
% axes('position',[0.55,0.55,0.3,0.3]);     %关键在这句！所画的小图
% plot(tl(1:100,1),cl(1:100,1),'g-','LineWidth',4.5); 
% % xlabel('t');
% % ylabel('y');
% %xlim([len-100,len]);%设置横坐标范围
% grid
% hold off
% 
for k=1:tt
    tr(1,k)=k;
end

figure(2)
subplot(3,1,1)
plot(tr(1,1:tt)', cen_C(1,1:tt)');figure(gcf);
hold on
ylabel('Statistics Times','FontSize',18);
xlabel(' ||L\Lambda ||_2','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectories of || L\Lambda||_2')
legend('||L\Lambda ||_2')
% %ylabel('The consensus error ||L\Lambda||_2^2 of five agents')
grid



subplot(3,1,2)

plot(tr(1,1:tt)',b_C(1,1:tt)');figure(gcf);

hold on 
ylabel('Statistics Times','FontSize',18);
xlabel('Balance gap ','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the power balance gap')
grid
z=legend('$$\sum  ||X-D||_2  $$');
set(z,'Interpreter','latex');



subplot(3,1,3)

plot(tr(1,1:tt)',op_C(1,1:tt)');figure(gcf);

 hold on
ylabel('Statistics Times','FontSize',18);
xlabel('Optimal criterion','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the optimal criterion')
grid
z= legend('$$|| X||_2 + || \dot{\Lambda}||_2 + ||\dot{Z} ||_2 $$');
set(z,'Interpreter','latex'); 

lent=len;

ll=100;
figure(3)
subplot(3,1,1)
tl=t(len-ll:len,1);
cl=censample(len-ll:len,1);
plot(tl(1:100,1),cl(1:100,1),'g-','LineWidth',3.5);
hold on
xlabel('Time','FontSize',18);
ylabel(' ||L\Lambda ||_2','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectories of || L\Lambda||_2')
legend('||L\Lambda ||_2')
% %ylabel('The consensus error ||L\Lambda||_2^2 of five agents')
grid
axes('position',[0.39,0.75,0.25,0.15]);     %关键在这句！所画的小图
 plot(t(1:len,1),censample(1:len,1),'g-','LineWidth',4.5)
% xlabel('t');
% ylabel('y');
%xlim([len-100,len]);%设置横坐标范围
grid
hold off

subplot(3,1,2)
tl=t(len-ll:len,1);
cl=bsample(len-ll:len,1);
plot(tl(1:100,1),cl(1:100,1),'b-','LineWidth',3.5); 
hold on 
xlabel('Time','FontSize',18);
ylabel('Balance gap ','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the power balance gap')
grid
z=legend('$$\sum  P_i^g- P_i^d $$');
set(z,'Interpreter','latex');


axes('position',[0.39,0.45,0.25,0.15]);     %关键在这句！所画的小图

plot(t(1:len,1),bsample(1:len,1),'b-','LineWidth',3.5)
grid
hold off

tt=1;
subplot(3,1,3)
tl=t(len-ll:len,1);
cl=opsample(len-ll:len,tt);
plot(tl(1:100,1),cl(1:100,1),'k-','LineWidth',3.5);
hold on
xlabel('Time','FontSize',18);
ylabel('Optimal criterion','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the optimal criterion')
grid
z= legend('$$|| \dot{P}^g||_2 + || \dot{\Lambda}||_2 + ||\dot{Z} ||_2 $$');
set(z,'Interpreter','latex'); 



axes('position',[0.39,0.15,0.25,0.15]);     %关键在这句！所画的小图
 
plot(t(1:len,1),opsample(1:len,tt),'k-','LineWidth',3.5)
grid
hold off

%% projection consensus for linear equation, the case of unique solution 
% The algorithm is over a random graph 
% yipeng 2020.02.28
% Problem: distributedly solving z=Hy
% The agent communication with random network G=(N,E), 
%  each agent has a few rows of H

%% Problem setting
% The linear equation problem z=Hy
clear all
clc
clf
% The number of agents N 
N=10; % This is  the number of the agent
% The dimension of the solution variable 
m=2;


% how many rows each agent holds
% The allocation vector, and the sum of rs should be r
rs=randi([1,5],1,N);
r=sum(rs); % This is the total number of row of H. H is a r*m matrix


%generate the data H for  The unique solution case
%Here we assume that rank (H)=m, and z\in span{h_1,h_2,...,\cdots,h_r}
%We first generate a random matrix r*r, which is a positive definite
%matrix, and  we take  m column, in this case, rank(H)=m     
Xd=diag(-5+10*rand(r,1)); 
flag1=0;
while flag1==0
Ut = orth(10*randn(r,r));
if rank(Ut)==r
    flag1=1;
end
end
    
HN = Ut' * Xd * Ut;
H=HN(:,1:m);
% generate the true y
ytrue1=-5+10*rand(m,1);

if rank(H)==m
    disp('Right matrix: H, unique solution z=Hy')
end
% The observation
z = H*ytrue1+ min(max(0.1*rand(r,1),-1),1);

% construct the projection matrix
DH=[];
zh=[];
 for i=1:N
       if i==1
           Ht=H(1:rs(1,i),:);
           zt=z(1:rs(1:i),:);
       else
       Ht=H( sum(rs(1:i-1))+1: sum( rs(1:i)), : );
       zt=z( sum(rs(1:i-1))+1: sum( rs(1:i)), : );
       end
       
     
       
       DH=blkdiag(DH, Ht'*Ht);
       zh=[zh; Ht'*zt];
 end
disp('The projection matrix is given')

%% generate NUMG random graph 
% random chose one during the algorithm 
% make sure that the average graph is check

NUMG=30;
% number of disconnected graph
NUM_DIS=0;
% store all the Laplacian matrix 
LT   = zeros(N,NUMG*N);
ADJT = zeros(N,NUMG*N);

% the average Laplacian
ADJS=zeros(N,N);
Lavg=zeros(N,N);

pmin=0.01;
pmax=0.05;

  flag=0; % For correctly generat NUMG graph 
  while flag==0
      NUM_DIS=0;
      for i=1:NUMG
       % the graph edge probability is random
       p  = rand*(pmax-pmin)+pmin;
       adj = randomGraph(N,p);
         if isConnected(adj)==0
          NUM_DIS=NUM_DIS+1;
        end  
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
 if egn(N-1,1)> 0+1000*eps
     disp('Average graph connected ! Right')
 else
     disp('Wrong communication graph!')
 end

disp('Create the graphs and weights')
 %% prepare for iteration
step=80000;

x1=zeros(N*m,step);
x2=zeros(N*m,step);
x3=zeros(N*m,step);
x4=zeros(N*m,step);

x1(:,1)=-5+10*rand(N*m,1);
x2(:,1)=x1(:,1);
x3(:,1)=x1(:,1);
x4(:,1)=x1(:,1);
% The algorithm
% x_{k+1}=DH£¨W \otimes I_m x_{k}£©+ zh
%initial value

% distance to optimal solution
dis1=zeros(step,1);
% consensus error
cen1=zeros(step,1);
% The error vector of the first agent
err11=zeros(m,step);

estimation1=zeros(1,step);

% distance to optimal solution
dis2=zeros(step,1);
% consensus error
cen2=zeros(step,1);
% The error vector of the first agent
err12=zeros(m,step);

estimation2=zeros(1,step);


% distance to optimal solution
dis3=zeros(step,1);
% consensus error
cen3=zeros(step,1);
% The error vector of the first agent
err13=zeros(m,step);

estimation3=zeros(1,step);


% distance to optimal solution
dis4=zeros(step,1);
% consensus error
cen4=zeros(step,1);
% The error vector of the first agent
err14=zeros(m,step);

estimation4=zeros(1,step);



ytrue=lsqr(H,z);
ytrue_m=zeros(m,step);
ytrue_m(:,1)=ytrue;


ytbar=kron(ones(N,1),ytrue);
ave=1/N*kron(ones(N,1)*ones(1,N),eye(m));


%% Generate a markov chain
p=0.1;
flag=0; % For correctly generat a connected graph 
  while flag==0
 adj = randomGraph(NUMG,p);     
 egn= graphSpectrum(adj);
 egn2=egn(NUMG-1,1); % whether the average graph is connected
 if egn(NUMG-1,1)> 0+eps
    flag=1;
 end
 end
L = laplacianMatrix(adj);
ht=1/(2*max(diag(L)));

P_Markov1 = eye(NUMG)-ht*L;
RW=1+rand(NUMG,NUMG);

P_Markov2=P_Markov1.*RW;
P_Markov=((1./sum(P_Markov2,2))*ones(1,NUMG)).*P_Markov2;


mc = dtmc(P_Markov);
figure(2);
graphplot(mc,'ColorEdges',true);

RGS_random = simulate(mc,step);
disp('Create the Markov chain')





beta=1;
tau1=0.5;
tau2=0.6;
tau3=0.65;
tau4=0.7;

alpha1=beta-tau1-0.01;
alpha2=beta-tau2-0.01;
alpha3=beta-tau3-0.01;
alpha4=beta-tau4-0.01;



gamma=1;  a0=1  ;
%% iteration
for k=1:step-1
    
% select one random graph from LT
     RGS   =  RGS_random(k,1);
     LR    =  LT(:,(RGS-1)*N+1:RGS*N);
     adjR  =  ADJT(:,(RGS-1)*N+1:RGS*N);    

gamma1=a0/(10+k)^(alpha1);
gamma2=a0/(10+k)^(alpha2);
gamma3=a0/(10+k)^(alpha3);
gamma4=a0/(10+k)^(alpha4);

delta=1/(1+k)^(beta);

% update
x1(:,k+1) = x1(:,k) - gamma1* kron(LR, eye(m))*x1(:,k) - delta*(DH*x1(:,k)- zh );
x2(:,k+1) = x2(:,k) - gamma2* kron(LR, eye(m))*x2(:,k) - delta*(DH*x2(:,k)- zh );
x3(:,k+1) = x3(:,k) - gamma3* kron(LR, eye(m))*x3(:,k) - delta*(DH*x3(:,k)- zh );
x4(:,k+1) = x4(:,k) - gamma4* kron(LR, eye(m))*x4(:,k) - delta*(DH*x4(:,k)- zh );

% performance
dis1(k,1)=(k+1)^(tau1)*norm(x1(:,k)-ytbar)^2;
cen1(k,1)=(k+1)^(tau1)*norm(x1(:,k)-ave*x1(:,k))^2;
err11(:,k)=x1(1:m,k)-ytrue;
ytrue_m(:,k)=ytrue;

dis2(k,1)=(k+1)^(tau2)*norm(x2(:,k)-ytbar)^2;
cen2(k,1)=(k+1)^(tau2)*norm(x2(:,k)-ave*x2(:,k))^2;
err12(:,k)=x2(1:m,k)-ytrue;

dis3(k,1)=(k+1)^(tau3)*norm(x3(:,k)-ytbar)^2;
cen3(k,1)=(k+1)^(tau3)*norm(x3(:,k)-ave*x3(:,k))^2;
err13(:,k)=x3(1:m,k)-ytrue;

dis4(k,1)=(k+1)^(tau4)*norm(x4(:,k)-ytbar)^2;
cen4(k,1)=(k+1)^(tau4)*norm(x4(:,k)-ave*x4(:,k))^2;
err14(:,k)=x4(1:m,k)-ytrue;

 end
disp('finish the iteration')
%% graph parts
linefont=4;
fontsize=24;

figure(3)
subplot(1,2,1)
semilogy(dis1(2:step-1,1)','r-', 'LineWidth', linefont)
hold on
semilogy(dis2(2:step-1,1)','k--', 'LineWidth', linefont)
semilogy(dis3(2:step-1,1)','b-.', 'LineWidth', linefont)
semilogy(dis4(2:step-1,1)','c*', 'LineWidth', linefont)
xlabel('Iteration','FontSize',24);
za=legend('$$\delta_1 = 0.49, \tau=0.5$$',...
'$$\delta_1 = 0.39, \tau=0.6$$',...
'$$\delta_1 = 0.34, \tau=0.65$$',...
'$$\delta_1 = 0.29, \tau=0.7$$',...
'FontSize',36);
set(za,'Interpreter','latex');
grid
z=ylabel('$$(k+1)^{\tau}\sum_{i=1}^{N}||x_{i}(k) -x^*||^2$$','FontSize',13);
set(z,'Interpreter','latex');
set(gca,'FontName','TimesNew Roman','FontSize',28); hold off


subplot(1,2,2)
semilogy(cen1(2:step-1,1)','r-', 'LineWidth', linefont)
hold on
semilogy(cen2(2:step-1,1)','k--', 'LineWidth', linefont)
semilogy(cen3(2:step-1,1)','b-.', 'LineWidth', linefont)
semilogy(cen4(2:step-1,1)','c*', 'LineWidth', linefont)

xlabel('Iteration','FontSize',24);
grid

za=legend('$$\delta_1 = 0.49, \tau=0.5$$',...
'$$\delta_1 = 0.39, \tau=0.6$$',...
'$$\delta_1 = 0.34, \tau=0.65$$',...
'$$\delta_1 = 0.29, \tau=0.7$$',...
'FontSize',36);
set(za,'Interpreter','latex');
z=ylabel('$$(k+1)^{\tau}\sum_{i=1}^{N} ||x_{i}(k) - \frac{1}{N}\sum_{i=1}^{N} x_{i}(k) ||^2$$','FontSize',13);
set(z,'Interpreter','latex');
set(gca,'FontName','TimesNew Roman','FontSize',28); hold off
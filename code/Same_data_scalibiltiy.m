% test for distributed resource allocation algorithm
% min \sum_{i=1}^N  f_i(x_i)
%   \sum x_i = \sum d_i
clear all
% number of the agents
N=1000;
% the probability of the random graph
pmax=0.01;
pmin=0.005;
%  Euler step-size and algorithm len
es=0.01;
len=1000;

tic
% test number

CA=zeros(N,1); % for a_i

CA=2.+ 4.*rand(N,1);

BA=zeros(N,1);

BA=1.+ 5.*rand(N,1);

% Generation the load demand vector

D=zeros(N,1);

D=normrnd(6,6,[N,1]);

TD=sum(D);

if TD>=0
    disp('Demand positive');
else
    disp('Demand negative');
end


tt=20;
% simulation results
 op_C=zeros(1,tt); % optimal norm 
 b_C=zeros(1,tt);  % power balance
 ds_R=zeros(1,tt); % the desity of the graph
 se_R=zeros(1,tt); % the second eigenvalue
opsample=zeros(len,tt);
for tnm=1:1:tt

% Generate the cost functions 
% f_i(x_i)= a_i x_i^2 +  b_i x_i 

% generation the random graph with Laplacian
  flag=0; % connected flag
  p=tnm*(pmax-pmin)/tt+pmin;
  while flag==0
      adj=randomGraph(N,p);
      L= laplacianMatrix(adj);
      if isConnected(adj)==1
          flag=1;
      else 
          flag=0;
      end            
  end
  L;
  adj;
  density = linkDensity(adj);
  ds_R(1,tnm)=density;
  % display the egnvalue
 egn= graphSpectrum(adj);
 %disp eignvalue
 disp('second eigen')
 egn2=egn(N-1,1)
 if egn(N-1,1)> 0+eps
     disp('Right')
 else
     disp('Wrong!')
 end
 se_R(1,tnm)=egn2;
 % algorithm part
 
 %inition
 x=zeros(N,len+1);
 lam=zeros(N,len+1);
 z=zeros(N,len+1);
 
 % algorithm part
 x(:,1)=D;
 for k=1:1:len
     %gradient
     gg=CA.*x(:,k)+BA;
     %iterations
     x(:,k+1)   =  x(:,k)   + es.* ( - gg + lam(:,k) ); 
     
     lam(:,k+1) = lam(:,k)  + es.*(-L*lam(:,k)- L*z(:,k)+ (D-x(:,k)));
     
     z(:,k+1)   =  z(:,k)   + es.* (L*lam(:,k));
     opsample(k,tnm)= norm( ( - gg + lam(:,k) ) )+ norm( -L*lam(:,k)- L*z(:,k)+ (D-x(:,k))   ) + norm((L*lam(:,k)));
 end
 
 % test the  optimal 
 gg=CA.*x(:,len+1)+BA;
 test1 = norm ( - gg + lam(:,len+1) );     
 test2 = norm (-L*lam(:,len+1)- L*z(:,len+1)+ (D-x(:,len+1)));
 test3 = norm (L*lam(:,len+1));
 
 disp('test optimal')
 test=test1+test2+test3;
 if abs(test)< 20
    disp('optimal Sucess')
 else
     disp('optimal fail')
 end
 op_C(1,tnm)=test;
 
    
 %test the equality 
 tg=sum(x(:,len+1));
 td=sum(D);
 disp('Balance')
 balance=abs(tg-td)/td;
 if balance <0.01
     disp('balance sucesss')
 else
     disp('balance fail')
 end
 b_C(1,tnm)=balance;
 
 
end
 

toc 

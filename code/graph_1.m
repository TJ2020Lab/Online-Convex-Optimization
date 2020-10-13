%% The graph part show the final result
clc
clf
for k=1:len
    t(k,1)=k;
end
tt=500;
for k=1:tt
    tr(1,k)=k;
end
bin=20;
linar=30;

figure(1)
subplot(3,1,1)
scen_C=sort(cen_C(1,1:tt));
hist(scen_C(1,1:tt-linar),bin);figure(gcf);
hold on
ylabel('Statistics Times','FontSize',18);
% z=xlabel('$$||\bar{L}\otimes I_m \Lambda ||_2$$');
% set(z,'Interpreter','latex');
xlabel('Consensus error','FontSize',18);
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectories of || L\Lambda||_2')
z=legend('$$||\bar{L}\otimes I_m\Lambda ||_2$$')
set(z,'Interpreter','latex');
% %ylabel('The consensus error ||L\Lambda||_2^2 of five agents')
grid
hold off

subplot(3,1,2)
sb_C=sort(b_C(1,1:tt));
hist(sb_C(1,1:tt-linar),bin);figure(gcf);

hold on 
ylabel('Statistics Times','FontSize',18);
xlabel('Resource Constraints ','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the power balance gap')
grid
z=legend('$$||1_n^T\otimes I_m (X-D)||_2  $$');
set(z,'Interpreter','latex');
hold off


subplot(3,1,3)
sop_C=sort(op_C(1,1:tt));
hist(sop_C(1,1:tt-linar),bin);figure(gcf);

 hold on
ylabel('Statistics Times','FontSize',18);
xlabel('Optimal criterion','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the optimal criterion')
grid
z= legend('$$ || P(S)||_2 $$');
set(z,'Interpreter','latex'); 

hold on
hold off


%%
lent=len;

ll=100;
figure(3)

subplot(3,1,1)

tl=t(len-ll:len,1);
cl=censample(len-ll:len,1);
plot(tl(1:100,1),cl(1:100,1),'g-','LineWidth',3.5);
hold on
xlabel('Recursive Times','FontSize',18);
z=ylabel('$$ ||\bar{L}\otimes I_m\Lambda ||_2$$')
set(z,'Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',20);

%title('The trajectories of || L\Lambda||_2')
z=legend('$$||\bar{L}\otimes I_m \Lambda ||_2$$')
set(z,'Interpreter','latex');
% %ylabel('The consensus error ||L\Lambda||_2^2 of five agents')
grid
axes('position',[0.39,0.75,0.25,0.15]);     %关键在这句！所画的小图
 plot(t(1:len,1),censample(1:len,1),'r-','LineWidth',4.5)
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
xlabel('Recursive Times','FontSize',18);
ylabel('Resource Constraint','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the power balance gap')
grid
z=legend('$$ ||1_n^T\otimes I_m (X-D) || $$');
set(z,'Interpreter','latex');


axes('position',[0.39,0.45,0.25,0.15]);     %关键在这句！所画的小图

plot(t(1:len,1),bsample(1:len,1),'r-','LineWidth',3.5)
grid
hold off


subplot(3,1,3)
tl=t(len-ll:len,1);
cl=opsample(len-ll:len,1);
plot(tl(1:100,1),cl(1:100,1),'k-','LineWidth',3.5);
hold on
xlabel('Recursive Times','FontSize',18);
ylabel('Optimal criterion','FontSize',24)
set(gca,'FontName','Times New Roman','FontSize',20);
%title('The trajectory of the optimal criterion')
grid
z= legend('$$ ||P(S) ||_2 $$');
set(z,'Interpreter','latex'); 



axes('position',[0.39,0.15,0.25,0.15]);     %关键在这句！所画的小图
 
plot(t(1:len,1),opsample(1:len,1),'r-','LineWidth',3.5)
grid
hold off

%% plot the graph
for i=1:N
      for j=1:N
          if abs(ADJS(i,j)) > 100*eps
              tempadj(i,j)=1;
          end
      end
end
figure(4)
subplot(2,2,4)
draw_circ_graph(tempadj)
hold on
xlabel('The Mean value graph','FontSize',18);
hold off

subplot(2,2,2)
adj=ADJT(:,3*N+1:4*N);
draw_circ_graph(adj)
hold on
xlabel('Sample Graph','FontSize',18);
hold off

subplot(2,2,3)
adj=ADJT(:,6*N+1:7*N);
draw_circ_graph(adj)
hold on
xlabel('Sample Graph','FontSize',18);
hold off

subplot(2,2,1)
adj=ADJT(:,10*N+1:11*N);
draw_circ_graph(adj)
hold on
xlabel('Sample Graph','FontSize',18);
hold off

% figure(2)
% subplot(3,1,1)
% plot(tr(1,1:tt)', cen_C(1,1:tt)');figure(gcf);
% hold on
% ylabel('Statistics Times','FontSize',18);
% xlabel(' ||L\Lambda ||_2','FontSize',24);
% set(gca,'FontName','Times New Roman','FontSize',20);
% %title('The trajectories of || L\Lambda||_2')
% legend('||L\Lambda ||_2')
% % %ylabel('The consensus error ||L\Lambda||_2^2 of five agents')
% grid
% hold off
% 
% subplot(3,1,2)
% 
% plot(tr(1,1:tt),b_C(1,1:tt));figure(gcf);
% 
% hold on 
% ylabel('Statistics Times','FontSize',18);
% xlabel('Balance gap ','FontSize',24)
% set(gca,'FontName','Times New Roman','FontSize',20);
% %title('The trajectory of the power balance gap')
% grid
% z=legend('$$\sum  ||X-D||_2  $$');
% set(z,'Interpreter','latex');
% hold off
% 
% 
% subplot(3,1,3)
% 
% plot(tr(1,1:tt),op_C(1,1:tt));figure(gcf);
% 
%  hold on
% ylabel('Statistics Times','FontSize',18);
% xlabel('Optimal criterion','FontSize',24)
% set(gca,'FontName','Times New Roman','FontSize',20);
% %title('The trajectory of the optimal criterion')
% grid
% z= legend('$$|| X||_2 + || \dot{\Lambda}||_2 + ||\dot{Z} ||_2 $$');
% set(z,'Interpreter','latex'); 
% 
% hold on
% hold off
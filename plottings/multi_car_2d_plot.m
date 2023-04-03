%% plotting (2D obstacle avoidance)
function multi_car_2d_plot(T,X,x0,xf,obstacles)
figure(1)
% subplot(2,2,1);
for ii=1:size(obstacles,1)
   rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on
end
scatter(xf(1),xf(2),80,'filled','hg'); hold on
scatter(x0(1),x0(2),50,'filled','or');
scatter(xf(5),xf(6),80,'filled','hb'); hold on
scatter(x0(5),x0(6),50,'filled','oy');

plot(X(1,:),X(2,:),'LineWidth',0.5); hold on;
plot(X(5,:),X(6,:),'LineWidth',0.5); hold on;
% title('x_1 vs x_2'); 
xlabel('$x$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
ylabel('$y$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% box on; axis square;
box on; axis equal;

% figure(2) % plot BaS
% plot(T,X(end,:),'-','LineWidth',0.5); hold on;
% % title('BaS vs time');
% ylabel('Barrier state $w_k$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% xlabel('Time ($s$)','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% legend('$(x_0,y_0)=(-3,-4)$','$(x_0,y_0)=(-4,-2)$','$(x_0,y_0)=(-2,-0.5)$','$(x_0,y_0)=(-4,2)$','Interpreter','latex');
% box on;

% figure(3)
% subplot(2,2,2)
% plot(T,X(1,:),'--',T,X(2,:),':'); title('x_1, x_2 vs t'); legend('x_1','x_2');
% subplot(2,2,3)
% plot(T(2:end),U(1,:),'--',T(2:end),U(2,:),':'); title('u_1, u_2 vs t'); legend('u_1','u_2');
% subplot(2,2,4)
% figure(4)
% plot(J); title('Cost vs iteration');
end
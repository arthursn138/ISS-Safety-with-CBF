%% plotting (2D obstacle avoidance)
function car_2d_plot(T,X,U,x0,xf,u_lims,obstacles)
% % % figure(1)
% % % subplot(2,2,1);
% % % for ii=1:size(obstacles,1)
% % %    rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on
% % % end
% % % scatter(xf(1),xf(2),80,'filled','hg'); hold on
% % % scatter(x0(1),x0(2),50,'filled','or');
% % % plot(X(1,:),X(2,:),'LineWidth',0.5); hold on;
% % % title('x_1 vs x_2'); 
% % % xlabel('$x$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% % % ylabel('$y$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% % % box on; axis square;
% % % box on; axis equal;
% % % 
% % % subplot(2,2,2)
% % % plot(T,X(1,:),'--',T,X(2,:),':'); title('x_1, x_2 vs t'); legend('x_1','x_2');
% % % 
% % % subplot(2,2,3)
% % % plot(T(2:end),U(1,:),'--',T(2:end),U(2,:),':'); title('u_1, u_2 vs t'); legend('u_1','u_2');
% % % % subplot(2,2,4)
% % % % figure(4)
% % % % plot(J); title('Cost vs iteration');

% % % figure % Plot BaS
% % % plot(T,X(end,:)); hold on; grid on;
% % % title('BaS vs time','Interpreter','latex');
% % % ylabel('Barrier state $w_k$','FontName','Times New Roman','Interpreter','latex');
% % % xlabel('Time ($s$)','FontName','Times New Roman','Interpreter','latex');
% % % % legend('$(x_0,y_0)=(-3,-4)$','$(x_0,y_0)=(-4,-2)$','$(x_0,y_0)=(-2,-0.5)$','$(x_0,y_0)=(-4,2)$','Interpreter','latex');
% % % box on;

% % % figure
% % % % subplot(2,2,1);
% % % plot(X(1,:),X(2,:),'LineWidth',0.5); hold on;
% % % title('x $(x_{1})$ $\times$ y $(x_{2})$', 'FontName', 'Times New Roman', 'Interpreter', 'latex'); 
% % % xlabel('$x$', 'FontName', 'Times New Roman', 'Interpreter', 'latex');
% % % ylabel('$y$', 'FontName', 'Times New Roman', 'Interpreter', 'latex');
% % % % box on; axis square;
% % % axis equal;
% % % 
% % % % subplot(2,2,2)
% % % figure
% % % plot(T, X(1,:));
% % % hold on; grid on;
% % % plot(T, X(2,:));
% % % title('States over time'); 
% % % legend('x $(x_{1})$', 'y $(x_{2})$', 'FontName', 'Times New Roman', 'Interpreter', 'latex');
% % % ylim([min(min(X(1,:),X(2,:))) 1.05*max(max(X(1,:),X(2,:)))]);
% % % 
% % % % subplot(2,2,3)
% % % figure
% % % plot(T(2:end), U(1,:));
% % % hold on; grid on;
% % % plot( T(2:end), U(2,:));
% % % title('Controls over time'); 
% % % legend('$u_{1}$', '$u_{2}$', 'FontName', 'Times New Roman', 'Interpreter', 'latex');
% % % ylim([(min(min(U(1,:),U(2,:))))-0.05 1.15*max(max(U(1,:),U(2,:)))]);
% % % 
% % % % subplot(2,2,4)
% % % % plot(J); title('Cost vs iteration');

%%%% ==== Yet another option (Arthur)
figure
subplot(2,2,1)

subplot(2,2,1);
for ii=1:size(obstacles,1)
   rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on
end
scatter(xf(1),xf(2),80,'filled','hg'); hold on
scatter(x0(1),x0(2),50,'filled','or');
plot(X(1,:),X(2,:),'LineWidth',1); hold on;
title('X-Y path (top view)','Interpreter','latex'); 
ylabel('$y$ (m)','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$x$ (m)','FontName','Times New Roman', 'Interpreter','latex');
box on; axis square; axis equal;

subplot(2,2,2)
plot(T,X(1,:),'LineWidth',1.5); hold on; grid on;
plot(T(end)+(T(end)-T(end-1)),xf(1),'*r','LineWidth',1.5)
% % title('X position vs time','Interpreter','latex');
% % ylabel('$x$ (m)','FontName','Times New Roman','Interpreter','latex');
% % xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');

plot(T,X(2,:),'LineWidth',1.5);
plot(T(end)+(T(end)-T(end-1)),xf(2),'*r','LineWidth',1.5)
plot(T,X(3,:),'LineWidth',1.5); hold on; grid on;
plot(T,xf(3)*ones(length(T)),'--r','LineWidth',1)
title('States over time','Interpreter','latex');
ylabel('$x$(m), $y$(m) and $\theta$(rad)','FontName','Times New Roman','Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
legend('$x$', '$x_{target}$', '$y$', '$y_{target}$', '$\theta$',...
    '$\theta_{target}$', 'FontName', 'Times New Roman','Interpreter','latex')
% box on;

subplot(2,2,3)
% % plot(T,X(3,:),'LineWidth',1.5); hold on; grid on;
% % plot(T,xf(3)*ones(length(T)),'--r','LineWidth',1.2)
% % title('Car heading vs time','Interpreter','latex');
% % ylabel('$\theta$ (rad)','FontName','Times New Roman','Interpreter','latex');
% % xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
plot(T,X(end,:)); hold on; grid on;
title('BaS vs time','Interpreter','latex');
ylabel('Barrier state $w_k$','Interpreter','latex');
xlabel('Time ($s$)','Interpreter','latex');
box on;

subplot(2,2,4)
% % plot(X(1,:),X(2,:),'LineWidth',1.5);hold on; grid on;
% % plot(xf(1),xf(2),'*r','LineWidth',1.5);
% % title('X-Y path (top view)','Interpreter','latex');
% % ylabel('$y$ (m)','FontName','Times New Roman', 'Interpreter','latex');
% % xlabel('$x$ (m)','FontName','Times New Roman', 'Interpreter','latex');
% % % box on;
% % axis equal
% % axis([-2 2 0 4])
% % 
% % figure
plot(T(1:end-1),U(1,:),'LineWidth',1.5); hold on; grid on;
plot(T(1:end-1),U(2,:),'LineWidth',1.5); 
title('Control vs time','FontName','Times New Roman','Interpreter','latex');
ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
legend('$u_{1}$', '$u_{2}$','FontName','Times New Roman','Interpreter','latex');
if ~isempty(u_lims)
    plot(T,u_lims(1,1)*ones(size(T)),'-.k','LineWidth',1.2);
    plot(T,u_lims(1,2)*ones(size(T)),'-.k','LineWidth',1.2);
    plot(T,u_lims(2,1)*ones(size(T)),':k','LineWidth',1.2);
    plot(T,u_lims(2,2)*ones(size(T)),':k','LineWidth',1.2);
    legend('$u_{1}$', '$u_{2}$','$u_{1_{min}}$', '$u_{1_{max}}$','$u_{2_{min}}$', '$u_{2_{max}}$',...
        'FontName','Times New Roman','Interpreter','latex');
end
if ~isempty(u_lims) && u_lims(1,1) == u_lims(2,1) && u_lims(1,2) == u_lims(2,2)
    ylim([1.1*u_lims(1,1) 1.1*u_lims(1,2)]);
end

% box on;

end
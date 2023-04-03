%% Cart-pole plotting
function cart_pole_plot(T,X,U,x0,xf,u_lims,cart_limit)
% % figure
% % % plot(T,X(1,:),'Color','#0072BD','LineWidth',1.5); hold on;
% % plot(T,X(1,:),'LineWidth',1.5); hold on;
% % title('Cart position vs time','Interpreter','latex');
% % ylabel('$x$ (m)','FontName','Times New Roman','Fontsize',12,'Interpreter','latex');
% % xlabel('Time (s)','FontName','Times New Roman','Fontsize',12,'Interpreter','latex');
% % % plot(T,cart_limit*ones(length(T),1),'--k',T,-cart_limit*ones(length(T),1),'--k','LineWidth',1.5);hold on;
% % % legend('Wall','DBaS-DDP','Unconstrained DDP','FontName','Times New Roman','Interpreter','latex','Fontsize',10);
% % box on;
% % 
% % figure
% % % plot(T,X(2,:),'Color','#0072BD','LineWidth',1.5); hold on;
% % plot(T,X(2,:),'LineWidth',1.5); hold on;
% % % plot(T,xf(2)*ones(length(T),1),'-g','LineWidth',1.5);
% % title('Pendulum angle vs time','Interpreter','latex');
% % ylabel('$\theta$ (rad)','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% % xlabel('Time (s)','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% % % legend('Targe angle','DBaS-DDP','Unconstrained DDP','FontName','Times New Roman','Fontsize',10);
% % box on;
% % 
% % figure
% % plot(T,X(3,:),'LineWidth',1.5);hold on;
% % title('Cart velocity vs time','Interpreter','latex');
% % ylabel('State x_3','FontName','Times New Roman');
% % xlabel('Time (s)','FontName','Times New Roman');
% % 
% % figure
% % plot(T,X(4,:),'LineWidth',1.5);
% % title('Pendulum velocity vs time','Interpreter','latex');
% % ylabel('$\dot{\theta}$ (rad/s)','FontName','Times New Roman', 'Interpreter','latex');
% % xlabel('Time (s)','FontName','Times New Roman');
% % box on;


% % Arthur's addition - all states in a single plot
zero = zeros(length(X));
pis = pi.*ones(length(X));

figure
subplot(3,2,1)
plot(T,X(1,:),'LineWidth',1.5); hold on; grid on;
plot(T,zero,'--r','LineWidth',1);
plot(T,cart_limit*ones(length(T),1),'k','LineWidth',3);
plot(T,-cart_limit*ones(length(T),1),'k','LineWidth',3);
title('Cart position vs time','Interpreter','latex');
ylabel('$x$ (m)','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
% legend('x', '$x_{target}$','Interpreter','latex');
box on;
subplot(3,2,2)
plot(T,X(2,:),'LineWidth',1.5); hold on; grid on;
plot(T,pis,'--r','LineWidth',1.5);
title('Pendulum angle vs time','Interpreter','latex');
ylabel('$\theta$ (rad)','FontName','Times New Roman','Fontsize',12,'Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
% box on;
subplot(3,2,3)
plot(T,X(3,:),'LineWidth',1.5); hold on; grid on;
plot(T,zero,'--r','LineWidth',1.5);
title('Cart velocity vs time','Interpreter','latex');
ylabel('$\dot{x}$ (m/s)','FontName','Times New Roman','Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
subplot(3,2,4)
plot(T,X(4,:),'LineWidth',1.5); hold on; grid on;
plot(T,zero,'--r','LineWidth',1.5);
title('Pendulum velocity vs time','Interpreter','latex');
ylabel('$\dot{\theta}$ (rad/s)','FontName','Times New Roman', 'Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman');
% box on;

% figure
% % % % % plot(T(1:end-1),U,'Color','#0072BD','LineWidth',1.5); hold on;
subplot(3,2,5)
plot(T(1:end-1),U,'LineWidth',1.5); hold on; grid on;
title('Control vs time','FontName','Times New Roman','Interpreter','latex');
ylabel('$u$ (N)','FontName','Times New Roman','Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
if ~isempty(u_lims)
    plot(T,u_lims(1,1)*ones(size(T)),'-.k','LineWidth',1.2);
    plot(T,u_lims(1,2)*ones(size(T)),'-.k','LineWidth',1.2);
    legend('$u$', '$u_{min}$', '$u_{max}$','FontName','Times New Roman','Interpreter','latex');
end
% plot(T,zero,'--r')
ylim([1.2*u_lims(1,1) 1.2*u_lims(1,2)]);
box on;

% figure
subplot(3,2,6)
plot(T,X(end,:),'LineWidth',1); grid on
title('BaS vs time','Interpreter','latex');
ylabel('Barrier state $w_k$','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
box on;

% % % % figure
% % % % subplot(3,2,1);
% % % % plot(T,X(1,:),'LineWidth',1.5); hold on;
% % % % % % % % plot(T,cart_limit.*ones(length(T)),'--k',T,-cart_limit.*ones(length(T)),'--k','LineWidth',1.5);
% % % % % % % % hold on;
% % % % % % % % % hold off;
% % % % title('Cart position vs time');
% % % % ylabel('State x_1','FontName','Times New Roman');
% % % % xlabel('Time (s)','FontName','Times New Roman');
% % % % subplot(3,2,2);
% % % % plot(T,X(2,:),'LineWidth',1.5); hold on;
% % % % plot(T,xf(2)*ones(length(T)),'LineWidth',1.5);
% % % % title('Pendulum angle vs time');
% % % % ylabel('State x_2','FontName','Times New Roman');
% % % % xlabel('Time (s)','FontName','Times New Roman');
% % % % subplot(3,2,3);
% % % % plot(T,X(3,:),'LineWidth',1.5);hold on;
% % % % title('Cart velocity vs time');
% % % % ylabel('State x_3','FontName','Times New Roman');
% % % % xlabel('Time (s)','FontName','Times New Roman');
% % % % subplot(3,2,4);
% % % % plot(T,X(4,:),'LineWidth',1.5);hold on;
% % % % title('Pendulum velocity vs time');
% % % % ylabel('State x_4','FontName','Times New Roman');
% % % % xlabel('Time (s)','FontName','Times New Roman');
% % % % subplot(3,2,5);
% % % % plot(T,X(end,:),'LineWidth',1.5);hold on;
% % % % title('BaS vs time');
% % % % ylabel('Barrier state w_k','FontName','Times New Roman');
% % % % xlabel('Time (s)','FontName','Times New Roman');
% % % % subplot(3,2,6);
% % % % plot(T(1:end-1),U,'LineWidth',1.5);hold on;
% % % % title('Control vs time');
% % % % ylabel('Input u','FontName','Times New Roman');
% % % % xlabel('Time (s)','FontName','Times New Roman');
end
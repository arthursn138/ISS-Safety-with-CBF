%% inv pend plotting
function inv_pend_plot(T,X,U,x0,xf,vel_limit,J)
figure(100);
plot(J); hold on

figure(1);
subplot(1,2,1);
plot(T,X(1,:),':','LineWidth',2); hold on; grid on;
xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
ylabel('$\theta$','FontName','Times New Roman','Interpreter','latex');
% axis square
subplot(1,2,2);
plot(T,X(2,:),':','LineWidth',2); hold on; grid on;
xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
ylabel('$\dot{\theta}$','FontName','Times New Roman','Interpreter','latex');
axis square
% set(gcf, 'Position',  [2000, 200, 1000, 500]);

figure(200); plot(U); hold on; 

end

%% 
% x1plot=[];
% for ii=1:length(X_sim(1,:))
%     if X_sim(1,ii)< -pi
% %         X_sim(1,ii) 
%     x1plot(ii)= X_sim(1,ii) + 2*pi;
%     while x1plot(ii) < -pi
%         x1plot(ii)= x1plot(ii) + 2*pi;
%     end
% %     x1plot(ii)
%     else if X_sim(1,ii)> pi
%         x1plot(ii)= X_sim(1,ii) - 2*pi;
%         else
%             x1plot(ii)= X_sim(1,ii);
%         end
%     end
% end
% figure()
% plot(T,x1plot)

% T=0:dt:dt*N-1*dt;
% figure(1);
% subplot(1,4,1);
% plot(T,X_sim(1,:),'LineWidth',1.5); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\theta$','FontName','Times New Roman','Interpreter','latex');
% legend('Penalty-DDP on Model','Penalty-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% % legend('DBaS-DDP on Model','DBaS-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% 
% % axis square
% subplot(1,4,2);
% plot(T,X_sim(2,:),'LineWidth',1.5); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\dot{\theta}$','FontName','Times New Roman','Interpreter','latex');
% % axis square
% set(gcf, 'Position',  [200, 200, 800, 300]);
% plot(T,-cart_limit*ones(length(T)),'k--','LineWidth',1); hold on;
% legend('Penalty-DDP on Model','Penalty-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% % legend('DBaS-DDP on Model','DBaS-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% % 
% legend('Robust DBaS-DDP on Model','Robust DBaS-DDP on True System','DBaS-DDP on Model','DBaS-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% 
% subplot(1,4,3);
% plot(0:dt:dt*N-2*dt,U_sim,'LineWidth',1.5); hold on
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% legend('Penalty-DDP on Model','Penalty-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% % legend('DBaS-DDP on Model','DBaS-DDP on True System','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% 
% 
% for ii=1:length(K_u)
% ku(ii,1:2)=K_u(1,1:2,ii);
% end
% subplot(1,4,4);
% plot(T,ku,'LineWidth',1.5); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$K_u(t)$ (feedback gain)','FontName','Times New Roman','Interpreter','latex');
% legend('$K_u^{x_1}$','$K_u^{x_2}$','Velocity Limit','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% title('Feedback gain of penalty method','FontName','Times New Roman','Interpreter','latex');
% figure(100);
% % subplot(1,2,1);
% plot(T,X_sim(3,:),'LineWidth',1.5); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\theta$','FontName','Times New Roman','Interpreter','latex');


% axis square
% subplot(1,2,2);
% plot(T,X_sim(2,:),'LineWidth',1.5); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\dot{\theta}$','FontName','Times New Roman','Interpreter','latex');
% % axis square
% set(gcf, 'Position',  [200, 400, 800, 300]);
% % plot(T,-cart_limit*ones(length(T)),'k--','LineWidth',1); hold on;
% figure(1)
% subplot(4,1,1);
% plot(T,X(1,:),'LineWidth',1.5); hold on
% % plot(T,xf(1)*ones(length(T)),'--','LineWidth',1);
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\theta$','FontName','Times New Roman','Interpreter','latex');
% subplot(4,1,2);
% plot(T,X(2,:),'LineWidth',1.5); hold on;
% plot(T,-cart_limit*ones(length(T)),'k--','LineWidth',1); hold on;
% legend('DBaS-MinMaxDDP on Model','DBaS-MinMaxDDP on True System','Velocity Limit','FontName','Times New Roman','FontSize',7.5,'Interpreter','latex')
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\dot{\theta}$','FontName','Times New Roman','Interpreter','latex');
% subplot(4,1,3);
% plot(0:dt:dt*N-2*dt,U,'LineWidth',1.5); hold on
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% subplot(4,1,4);
% plot(T,X(3,:),'LineWidth',1.5); hold on;
% % plot(T,xf(2)*ones(length(T)),'--','LineWidth',1);
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('Barrier State','FontName','Times New Roman','Interpreter','latex');

% figure(2)
% plot(J); hold on;
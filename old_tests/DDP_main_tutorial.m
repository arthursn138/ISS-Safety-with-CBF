%% Safety Embedded Differential Dynamic Programming using Discrete Barrier States (DBaS)
% Hassan Almubarak - ACDS Lab @ Georgia Tech
% halmubarak@gatech.edu
% Last Update May/29/2021

% To run this main file (SafeDDP_main),
% 1. call the system's dynamics
% 2. generate the obstacle course and define the safe set function (h)
% 3. call DBaS_dyn to generate the DBaS dynamics
% 4. call Safety_Embedding_dynamics to augment the DBaS to the system's
% dynamics
% 5. Define DDP and optimization paramters and run vanilla ddp

% NOTE: to avoid pentrating the obstacles in some cases due to the discrete
% formulation, use disc_ddp_alg_penalty which penalizes the interior of the
% unsafe regions as done in the penalty methods. disc_ddp_alg_penalty takes
% h as an extra input to penalize the interiors of the obstacles
% This gives advantages of good planning from safety embedded ddp and the
% advantage of penalty methods. 

% clear workspace and close figures
clear
% close all
% clc

%% add paths
genral_path = '../hassan_safety_embedded_ddp/';
addpath(genpath(genral_path))  
%% initialize system parameters
dt=0.01; %sampling
N=600; % horizon
T=0:dt:dt*N-1*dt;

% [f,fx,fu,x,u]=Rocket_landing_dynamics4(dt,1);
% [f,fx,fu,x,u]=PlanDInt_dynamics(dt,1);
% [f,fx,fu,x,u]=Diff_wheel_dynamics(dt,1);
% [f,fx,fu,x,u]=inv_pend_dynamics(dt,1);
% [f,fx,fu,x,u] = quadrotor_dynamics(dt,1);
[f,fx,fu,x,u]=cart_pole_dynamics(dt,1);
n=length(x); %state dim 
m=length(u); %input dim

% specify initial state for the original system
x0=[0;0;0;0];
% specify final desired state for the original system
xf=[0;pi;0;0];

f_dyn.f=f; f_dyn.fx=fx; f_dyn.fu=fu;

sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
%% Quadratic costs (running cost and terminal cost)
Q=0*eye(n); 

R=.1*eye(m);
S=100*eye(n);
S(2,2)=10000;

run_cost=@(x,u,deriv_bool) run_quad_cost(x,u,Q,R,xf,deriv_bool);
term_cost=@(x,deriv_bool) terminal_quad_cost(x,xf,S,deriv_bool);
%% nominal input and state
ubar = 0.01*ones(m, N-1); % nominal control
xbar=[]; xbar(:,1) = x0;          % initial state
for k=1:N-1
    xbar(:,k + 1) = f_dyn.f(xbar(:, k), ubar(:, k));
end
%% optimization parameters (by Yuichiro Aoyama)
iter= 20;                % number of max iterations
toler= 1e-3;             % cost change 1e-3
% for regularization part by Yuichiro Aoyama
lambda = 1;              % initial value for lambda for regularization
dlambda= 1;              % initial value for dlambda for regularization
lambdaFactor = 1.6;      % lambda scaling factor
lambdaMax = 1e10;        % lambda maximum value
lambdaMin = 1e-6;        % below this value lambda = 0
opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda', dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax, 'lambdaMin', lambdaMin);
%% Discrete DDP
[X,U,J,lambda,dlambda,alpha,k_u,K_u] = disc_ddp_alg(f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par);
%% cartpole plottings
figure(1)
subplot(4,1,1)
plot(T,X(1,:),'LineWidth',1.5); hold on;
ylabel('Cart position $x$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time [$s$]','FontName','Times New Roman','Interpreter','latex');
subplot(4,1,2)
plot(T,X(2,:),'LineWidth',1.5); hold on;
plot(T,3.14*ones(N,1),'g--','LineWidth',1);
ylabel('Pole angle $\theta$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time [$s$]','FontName','Times New Roman','Interpreter','latex');
subplot(4,1,3)
plot(T,X(3,:),'LineWidth',1.5); hold on;
ylabel('Cart velocity $\dot{x}$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time [$s$]','FontName','Times New Roman','Interpreter','latex');
subplot(4,1,4)
plot(T,X(4,:),'LineWidth',1.5); hold on;
ylabel('Pole angular velocity $\dot{\theta}$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time [$s$]','FontName','Times New Roman','Interpreter','latex');

figure(2)
plot(0:dt:dt*N-2*dt,U,'LineWidth',1.5); hold on;
ylabel('Controller $u$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time [$s$]','FontName','Times New Roman','Interpreter','latex');

cart_position = X(1,:);
theta = X(2,:);
x_pole=[]; y_pole=[]; L =2;
for ii=1:N
x_pole(ii)=cart_position(ii)+L*cos(theta(ii)-pi/2);
y_pole(ii)=L*sin(theta(ii)-pi/2);
end

fig =figure(100); 
axis tight manual % this ensures that getframe() returns a consistent size
for idx = 1:10:N
    plot(cart_position(idx),0,'s',[cart_position(idx),x_pole(idx)],[0;y_pole(idx)],'-',x_pole(1:idx),y_pole(1:idx),'--',x_pole(idx),y_pole(idx),'o','LineWidth',1.5); 
    axis([-5 5 -5 5]);
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end


%% inv pend plotting
% T=0:dt:dt*N-1*dt;
% figure(3);
% subplot(1,2,1);
% plot(T,X_sim(1,:),':','LineWidth',2); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\theta$','FontName','Times New Roman','Interpreter','latex');
% % axis square
% subplot(1,2,2);
% plot(T,X_sim(2,:),':','LineWidth',2); hold on; grid on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$\dot{\theta}$','FontName','Times New Roman','Interpreter','latex');
% axis square
% set(gcf, 'Position',  [2000, 200, 1000, 500]);
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
%%
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

%% for quad - many sims 
% XX1=[]; XX2=[]; XX3=[]; BB=[];
% T=0:dt:dt*N-1*dt;
% for kk=1:500
%     [X_sim,U_sim]=simulate_ddp(f_dyn,k_u,K_u,x0,N,X,U);
% % plotting (3D obstacle avoidance)
% % figure(10000)
% % plot3(X(10,:),X(11,:),X(12,:),'b-',x0(10),x0(11),x0(12),'bo',xf(10),xf(11),xf(12),'gx','LineWidth',1); hold on;
% % plot3(X(10,:),X(11,:),X(12,:),'LineWidth',1); hold on;
% % scatter3(xf(10),xf(11),xf(12),50,'filled','og','LineWidth',1); hold on
% % scatter3(x0(10),x0(11),x0(12),50,'filled','or','LineWidth',1);
% % plot3(X_sim(10,:),X_sim(11,:),X_sim(12,:),'-','Color','#D95319','LineWidth',.25); hold on;
% % plot3(X_sim(10,:),X_sim(11,:),X_sim(12,:),'-','LineWidth',.25); hold on;
% % grid on;
% % 
% % % 
% %     num_obs=size(obs_info,1);
% %     [sphere_x,sphere_y,sphere_z]=sphere;
% %     x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
% %     r = obs_info(:,4);
% %     figure(10000)
% %     for ii=1:num_obs
% %     spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
% %     hold on;
% %     spheres.EdgeColor = 'k';
% %     spheres.FaceColor = '#A2142F';
% %     spheres.LineStyle = ':';
% %     spheres.FaceLighting = 'flat';
% %     spheres.FaceAlpha= 1;
% %     end
% %     axis equal;   
% % xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
% % ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
% % zlabel('$z$','FontName','Times New Roman','Interpreter','latex');
% % %     axis([-6 17 -6 6 -6 6])
% % box on; axis square;
% % box on; axis equal;
% 
% % figure(50)
% % subplot(3,1,1);
% % plot(0:dt:dt*N-2*dt,U_sim,'LineWidth',1.5); hold on
% % xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% % ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% % subplot(3,1,2);
% % plot(0:dt:dt*N-2*dt,V_sim,'LineWidth',1.5); hold on
% % xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% % ylabel('$v$','FontName','Times New Roman','Interpreter','latex');
% % subplot(3,1,3);
% % plot(T,X_sim(13,:),'LineWidth',1.5); hold on;
% % xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% % ylabel('Barrier State $w$','FontName','Times New Roman','Interpreter','latex');
% XX1(kk,:)=X_sim(10,:);
% XX2(kk,:)=X_sim(11,:);
% XX3(kk,:)=X_sim(12,:);
% BB(kk,:)=X_sim(13,:);
% end
% 
% %% Confidence region approx
% figure(1000);
% % for ii=1:size(XX1,1)
% %     plot3(XX1(ii,:),XX2(ii,:),XX3(ii,:),'-','Color','#D95319','LineWidth',.25); hold on;
% % end
%     num_obs=size(obs_info,1);
%     x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
%     r = obs_info(:,4); [sphere_x,sphere_y,sphere_z]=sphere;
%     for ii=1:num_obs
%     spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
%     hold on;
%     spheres.EdgeColor = 'k';
%     spheres.FaceColor = '#A2142F';
%     spheres.LineStyle = ':';
%     spheres.FaceLighting = 'flat';
%     spheres.FaceAlpha= 1;
%     end
%     axis equal;   
%     
% facecolor='#D95319'; alpha=0.01; edgecolor = 'none';
% linestyle = 'none';
% conf_perc=99/100; spc =1;
% zstar=chi2inv(conf_perc,3);  % 3 is for 3D 
% fig_num=1000;
% plot_conf_reg(fig_num,XX1,XX2,XX3,N,zstar,spc,facecolor,alpha,edgecolor,linestyle)
% 
% plot3(X(10,:),X(11,:),X(12,:),'-','Color','#D95319','LineWidth',1.5); hold on;
% plot3(x0(10),x0(11),x0(12),'ro','LineWidth',1)
% plot3(xf(10),xf(11),xf(12),'gx','LineWidth',2); hold on;
% 
% % title('Safety Embedded Robust DDP','FontName','Times New Roman','Interpreter','latex');
% xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
% ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
% zlabel('$z$','FontName','Times New Roman','Interpreter','latex');
% % legend('Start','Target','Undisturbed Trajectory','$95\%$ Confidence Region','Unsafe Region','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
% grid on;

% for ii=1:size(XX1,1)
% plot3(XX1(ii,:),XX2(ii,:),XX3(ii,:),'-','Color',facecolor,'LineWidth',.25); hold on;
% end
%%
% Confidence interval for barrier state:
% BaS_mean = mean(BB);
% BaS_std = std(BB);
% zstar=sqrt(chi2inv(conf_perc,1));
% % upper bound:
% BaS_upper=BaS_mean + zstar*BaS_std;
% BaS_lower=BaS_mean - zstar*BaS_std;
% figure(51)
% plot(T,X(13,:),'-','Color','#D95319','LineWidth',1.5); hold on;
% % plot(T,BaS_mean,T,BaS_upper,T,BaS_lower); hold on;
% conf_int = patch([T fliplr(T)],[min([BaS_upper;BaS_lower]) fliplr(max([BaS_upper;BaS_lower]))],'r')
% conf_int.EdgeColor = edgecolor;
% conf_int.LineStyle = linestyle;
% conf_int.FaceColor = facecolor;
% conf_int.FaceAlpha = 0.5;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('Barrier State','FontName','Times New Roman','Interpreter','latex');
% legend('Undisturbed Trajectory','$95\%$ Confidence Region','Unsafe Region','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)


% figure(51);
% for kk=1:size(BB,1)
%     plot(T,BB(kk,:)); hold on;
% end





%% plotting (2D obstacle avoidance)
% T=0:dt:dt*N-1*dt;
% figure(1000)
% % subplot(2,2,1);
% % for ii=1:size(circ,1)
% %    rectangle('Position',circ(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on
% % end
% % scatter(xf(1),xf(2),80,'filled','hg'); hold on
% scatter(x0(1),x0(2),50,'filled','or');
% plot(X(1,:),X(2,:),'LineWidth',1.5); hold on;
% % title('x_1 vs x_2'); 
% xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
% ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
% % box on; axis square;
% box on; axis equal;
% % subplot(2,2,2)
% % plot(T,X(1,:),'--',T,X(2,:),':'); title('x_1, x_2 vs t'); legend('x_1','x_2');
% % subplot(2,2,3)
% % plot(T(2:end),U(1,:),'--',T(2:end),U(2,:),':'); title('u_1, u_2 vs t'); legend('u_1','u_2');
% % subplot(2,2,4)
% % plot(J); title('Cost vs iteration');
% % figure(2)
% % plot(T,X(n,:),'-','LineWidth',1.5); hold on;
% % % title('BaS vs time');
% % ylabel('Barrier state $w_k$','FontName','Times New Roman','Interpreter','latex');
% % xlabel('Time ($s$)','FontName','Times New Roman','Interpreter','latex');
% % % legend('$(x_0,y_0)=(-3,-4)$','$(x_0,y_0)=(-4,-2)$','$(x_0,y_0)=(-2,-0.5)$','$(x_0,y_0)=(-4,2)$','Interpreter','latex');
% % box on;
%% rocket plotting
% T=0:dt:dt*(N-1);
% figure(1)
% subplot(1,2,1)
% plot(x0(2),x0(1),'or');
% hold on
% plot(xf(2),xf(1),'og');
% plot(X(2,:),X(1,:));
% figure(1)
% subplot(1,2,2);
% plot(T(1:end-1),U)
% subplot(1,3,3);
% plot(T,X(end,:))

% figure(2)
% plot3(x0(2),x0(3),x0(1),'or'); hold on;
% plot3(xf(2),xf(3),xf(1),'og');
% plot3(X(2,:),X(3,:),X(1,:));
% axis([-3 3 -3 3 0 4]); grid on
%%
% angle_gamma=[];
% angle_theta=[];
% Tnorm=[];
% angle_delta=[];
% w=[];
% for ii=1:length(X)-1
%     
%     w(ii)= norm([X(11,:);X(12,:);X(13,:)]);
% %     const2(ii) = X(1,ii) - tand(gamma_gs) * norm([X(2,ii);X(3,ii)]);
%     angle_gamma(ii) = atand(X(1,ii)/norm([X(2,ii);X(3,ii)]));
% 
% %     const3(ii) = 1 - 2 * (X(9,ii)^2 + X(10,ii)^2) - cosd(theta_max);
%     angle_theta(ii) = acosd(1 - 2 * (X(9,ii)^2 + X(10,ii)^2));
%     
%     Tnorm(ii) = norm([X(14,ii);X(15,ii);X(16,ii)]);
%     
%     angle_delta(ii) = acosd(X(14,ii)/norm([X(14,ii);X(15,ii);X(16,ii)]));
% end
% w(ii+1)=w(ii);
% angle_gamma(ii+1)=angle_gamma(ii);
% angle_theta(ii+1)=angle_theta(ii);
% Tnorm(ii+1)=Tnorm(ii);
% angle_delta(ii+1)=angle_delta(ii);
% 
% 
% figure(5)
% subplot(5,1,1)
% plot(T,w,T,w_max*ones(length(T))); ylabel('w');
% subplot(5,1,2)
% plot(T,angle_gamma,T,gamma_gs*ones(length(T))); ylabel('\gamma');
% subplot(5,1,3)
% plot(T,angle_theta,T,theta_max*ones(length(T))); ylabel('\theta');
% subplot(5,1,4)
% plot(T,Tnorm,T,Tmax*ones(length(T))); ylabel('T');
% subplot(5,1,5)
% plot(T,angle_delta,T,delta_max*ones(length(T))); ylabel('\delta');


% figure()
% plot(T,X(11,:), T, X(12,:), T, X(13,:))

X(:,end)
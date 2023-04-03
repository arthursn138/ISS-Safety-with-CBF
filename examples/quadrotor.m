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
clear all
pause(1);
% close all
% clc

%% add paths
% cost_path = '../costs/';
% addpath(genpath(cost_path))                               

% ddp_path = '../ddp_algorithm/';
% addpath(ddp_path)     

% sys_path = '../Systems Dynamics/';
% addpath(sys_path)     

parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory)); 
%% initialize system parameters
dt=0.01; %sampling
N=401; % quad horizon
% N=101; % horizon
T=dt*N;

[f,fx,fu,x,u]=quadrotor_dynamics(dt,1);

n=length(x); %state dim (total, with Bas)
m=length(u); %input dim

% specify initial state for the original system
x0=zeros(n,1);
x0(10)=10; x0(11)=0; x0(12)=-1;
% specify final desired state for the original system
xf=zeros(n,1);
xf(10)=-5; xf(11)=-3; xf(12)=2;
if length(xf)~=n || length(x0)~=n
    error('wrong dimention of boundary conditions');
end
f_dyn.f=f; f_dyn.fx=fx; f_dyn.fu=fu;


barrier_state_button = 1;
penalty_button = 0;
ddp_2nd_order =0;

[h,obs_info]=obs_manual_3D;
if barrier_state_button == 1       
% generate the DBaS given h and f
    [f_w,fx_w,fu_w,w0,wf,fxx_w,fxu_w,fuu_w]=DBaS_dyn(x,u,x0,xf,f_dyn,h,ddp_2nd_order);
% augment the bas
    [fbar,fbarx,fbaru,fbarxx,fbarxu,fbaruu,xbar0,xbarf,nbar]=Safety_Embedding_dynamics(x,u,x0,xf,f_dyn,w0,wf,f_w,fx_w,fu_w,fxx_w,fxu_w,fuu_w,1,ddp_2nd_order);
% redfine the variables to use vanilla DDP
    f_dyn.f=fbar; f_dyn.fx=fbarx; f_dyn.fu=fbaru;
    x0=xbar0; xf=xbarf;n=nbar;
end

sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
%% Quadratic costs (running cost and terminal cost)
% quad :
Q=0*eye(n); 
Q_bf =0.2;
R=1*1e-4*eye(m); 

S=1*eye(n);
S(10:12,10:12)=2*eye(3);
S_bf=0.8;

run_cost=@(x,u,deriv_bool) run_quad_cost(x,u,Q,R,xf,deriv_bool);

if penalty_button == 1
    beta = 0;
        for ii=1:length(h)
            beta = beta + 1/h{ii}(x);
        end
    beta = matlabFunction(beta,'Vars',{x});  
    beta_x= jacobian(beta(x),x); beta_x=matlabFunction(beta_x,'Vars',{x});
    beta_xx= jacobian(beta_x(x),x);beta_xx=matlabFunction(beta_xx,'Vars',{x});
    run_cost=@(x,u,deriv_bool) run_quad_cost_penalty(x,u,Q,R,xf,beta,beta_x,beta_xx,Q_bf,deriv_bool);
end
if barrier_state_button == 1    
    Q(n,n)=Q_bf;
    S(n,n)=S_bf;
    run_cost=@(x,u,deriv_bool) run_quad_cost(x,u,Q,R,xf,deriv_bool);
end

term_cost=@(x,deriv_bool) terminal_quad_cost(x,xf,S,deriv_bool);
%% nominal input and state
ubar = 0.0*ones(m, N-1); % nominal control
% ubar(1,:)=9.81*ones(1,N-1);
xbar=[]; xbar(:,1) = x0;          % initial state
for k=1:N-1
    xbar(:,k + 1) = f_dyn.f(xbar(:, k), ubar(:, k));
end
%% optimization parameters
iter= 500;                % number of max iterations
toler= 1e-3;             % cost change 1e-3
% for regularization part by Yuichiro Aoyama
lambda = 0;              % initial value for lambda for regularization
dlambda= 0;              % initial value for dlambda for regularization
lambdaFactor = 1e-5;      % lambda scaling factor
lambdaMax = 1e10;        % lambda maximum value
lambdaMin = 1e-6;        % below this value lambda = 0
opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda', dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax, 'lambdaMin', lambdaMin);
%% Discrete DDP
% [X,U,J,lambda,dlambda,alpha,k_u,K_u,ii,iter_succ,L] =disc_ddp_alg(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par);
[X,U,J,lambda,dlambda,alpha,k_u,K_u,ii,iter_succ] =disc_ddp_alg_penalty(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par,h);
%% quadrotor plots
plot3(x0(10),x0(11),x0(12),'ro','LineWidth',1); hold on;
plot3(xf(10),xf(11),xf(12),'gx','LineWidth',2); hold on;
plot3(X(10,:),X(11,:),X(12,:),'-','Color',facecolor,'LineWidth',1.5); hold on;
% plot_conf_reg(fig_num,XX1,XX2,XX3,N,zstar,spc,facecolor,0.3,edgecolor,linestyle);
num_obs=size(obs_info,1);
x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
r = obs_info(:,4); [sphere_x,sphere_y,sphere_z]=sphere;
    for ii=1:num_obs
    spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
    hold on;
    spheres.EdgeColor = 'k';
    spheres.FaceColor = '#A2142F';
    spheres.LineStyle = ':';
    spheres.FaceLighting = 'flat';
    spheres.FaceAlpha= 1;
    end
axis equal; 
%% Simulate with no V player and possible change in dynamics
monte_carol_num = 100;
XX1=[]; XX2=[]; XX3=[]; BB=[];
T=0:dt:dt*N-1*dt;
for kk=1:monte_carol_num
    [X_sim,U_sim]=simulate_quadrotor(f_dyn,k_u,K_u,x0,N,X,U,dt);
% plotting (3D obstacle avoidance)
% figure(10)
% plot3(X(10,:),X(11,:),X(12,:),'b-',x0(10),x0(11),x0(12),'bo',xf(10),xf(11),xf(12),'gx','LineWidth',1); hold on;
% plot3(X(10,:),X(11,:),X(12,:),'LineWidth',1); hold on;
% scatter3(xf(10),xf(11),xf(12),50,'filled','og','LineWidth',1); hold on
% scatter3(x0(10),x0(11),x0(12),50,'filled','or','LineWidth',1);
% plot3(X_sim(10,:),X_sim(11,:),X_sim(12,:),'-','Color','#0072BD','LineWidth',.15); hold on;
% plot3(X_sim(10,:),X_sim(11,:),X_sim(12,:),'-','LineWidth',.25); hold on;
% grid on;
% 
% % 
%     num_obs=size(obs_info,1);
%     [sphere_x,sphere_y,sphere_z]=sphere;
%     x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
%     r = obs_info(:,4);
% %     figure(10000)
%     for ii=1:num_obs
%     spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
%     hold on;
%     spheres.EdgeColor = 'k';
%     spheres.FaceColor = '#A2142F';
%     spheres.LineStyle = ':';
%     spheres.FaceLighting = 'flat';
%     spheres.FaceAlpha= 1;
%     end
%     axis equal;   grid on;
% xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
% ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
% zlabel('$z$','FontName','Times New Roman','Interpreter','latex');
% %     axis([-6 17 -6 6 -6 6])
% box on; axis square;
% box on; axis equal;

% figure(50)
% subplot(3,1,1);
% plot(0:dt:dt*N-2*dt,U_sim,'LineWidth',1.5); hold on
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% subplot(3,1,2);
% plot(0:dt:dt*N-2*dt,V_sim,'LineWidth',1.5); hold on
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$v$','FontName','Times New Roman','Interpreter','latex');
% subplot(3,1,3);
% plot(T,X_sim(13,:),'LineWidth',1.5); hold on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('Barrier State $w$','FontName','Times New Roman','Interpreter','latex');
XX1(kk,:)=X_sim(10,:);
XX2(kk,:)=X_sim(11,:);
XX3(kk,:)=X_sim(12,:);
BB(kk,:)=X_sim(13,:);
end
% blue: #0072BD, purple: #7E2F8E
% Confidence region approx
facecolor='#D95319'; % orange
% facecolor='#0072BD'; % blue
% facecolor='#EDB120'; % yellow
% facecolor='#7E2F8E'; % purple
% facecolor='#77AC30'; % green
% facecolor='#4DBEEE'; % sky blue
% facecolor='#A2142F'; % brick red

% facecolor='#0072BD'; % blue
conf_perc=95/100; spc =1;
alpha=0.01; edgecolor = 'none';
linestyle = 'none';
zstar=chi2inv(conf_perc,3);  % 3 is for 3D 
fig_num=6652370;
figure(fig_num)
plot3(x0(10),x0(11),x0(12),'ro','LineWidth',1); hold on;
plot3(xf(10),xf(11),xf(12),'gx','LineWidth',2); hold on;
plot3(X(10,:),X(11,:),X(12,:),'-','Color',facecolor,'LineWidth',1.5); hold on;
% plot_conf_reg(fig_num,XX1,XX2,XX3,N,zstar,spc,facecolor,0.3,edgecolor,linestyle);
num_obs=size(obs_info,1);
x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
r = obs_info(:,4); [sphere_x,sphere_y,sphere_z]=sphere;
    for ii=1:num_obs
    spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
    hold on;
    spheres.EdgeColor = 'k';
    spheres.FaceColor = '#A2142F';
    spheres.LineStyle = ':';
    spheres.FaceLighting = 'flat';
    spheres.FaceAlpha= 1;
    end
axis equal;   
plot_confidence_region(fig_num,XX1,XX2,XX3,N,zstar,spc,facecolor,alpha,edgecolor,linestyle);

% title('Safety Embedded Robust DDP','FontName','Times New Roman','Interpreter','latex');
xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
zlabel('$z$','FontName','Times New Roman','Interpreter','latex');
% legend('Start','Target','Undisturbed Trajectory','$95\%$ Confidence Region','Unsafe Region','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
grid on;

% for ii=1:size(XX1,1)
% plot3(XX1(ii,:),XX2(ii,:),XX3(ii,:),'-','Color',facecolor,'LineWidth',.25); hold on;
% end
%%
% Confidence interval for barrier state:
BaS_mean = mean(BB);
BaS_std = std(BB);
zstar=sqrt(chi2inv(conf_perc,1));
% upper bound:
BaS_upper=BaS_mean + zstar*BaS_std;
BaS_lower=BaS_mean - zstar*BaS_std;
figure(51)
plot(T,X(13,:),'-','Color','#0072BD','LineWidth',1.5); hold on;
% plot(T,BaS_mean,T,BaS_upper,T,BaS_lower); hold on;
conf_int = patch([T fliplr(T)],[min([BaS_upper;BaS_lower]) fliplr(max([BaS_upper;BaS_lower]))],'r')
conf_int.EdgeColor = edgecolor;
conf_int.LineStyle = linestyle;
conf_int.FaceColor = facecolor;
conf_int.FaceAlpha = 0.5;
xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
ylabel('Barrier State','FontName','Times New Roman','Interpreter','latex');
legend('Undisturbed Trajectory','$95\%$ Confidence Region','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)

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
% cost_path = '../costs/';
% addpath(genpath(cost_path))                               

% ddp_path = '../ddp_algorithm/';
% addpath(ddp_path)     

% sys_path = '../Systems Dynamics/';
% addpath(sys_path)     

genral_path = '../hassan_safety_embedded_ddp/';
addpath(genpath(genral_path))  
%% initialize system parameters
dt=0.02; %sampling
N=500; % horizon
T=dt*N
T=0:dt:dt*N-1*dt;

% [f,fx,fu,x,u]=Rocket_landing_dynamics4(dt,1);
% m_wet=2;r_I_i=[2;0;2]; v_I_i=[-1;0.2;0]; q_B_I_i=[1;0;0;0]; w_B_i=[0;0;0];
% m_f=2;r_I_f=0e-3*ones(3,1); v_I_f=-1e-1*[1;0;0]; q_B_I_f=[1;0;0;0]; w_B_f=[0;0;0];
% [f,fx,fu,x,u]=PlanDInt_dynamics(dt,1);
[f,fx,fu,fxx,fxu,fuu,x,u]=Diff_wheel_dynamics(dt,1);
% [f,fx,fu,fxx,fxu,fuu,x,u]=cart_pole_dynamics(dt,1);
% [f,fx,fu,x,u]=inv_pend_dynamics(dt,1);
% [f,fx,fu,x,u] = quadrotor_dynamics(dt,1);
% [f,fx,fu,x,u] = toy_nonaffine_dynamics(dt,1);

n=length(x); %state dim (total, with Bas)
m=length(u); %input dim

% specify initial state for the original system
% x0 = zeros(n,1)+rand(3,1)/5;
x0 = zeros(n,1);
% x0=[-3.5;0;0];
% x0(10)=10; x0(11)=0; x0(12)=-1;
% x0=[pi;0];
% specify final desired state for the original system
% xf = zeros(n,1); % make sure to put the right dimention
% xf(10)=-5; xf(11)=-3; xf(12)=2;
% xf=[0;0];
xf=[3.5;3;0];
% xf=[0;pi;0;0];
if length(xf)~=n || length(x0)~=n
    error('wrong dimention of boundary conditions');
end
barrier_state_button = 1;
penalty_button = 0;
ddp_2nd_order = 0;

f_dyn.f=f; f_dyn.fx=fx; f_dyn.fu=fu; 
if ddp_2nd_order
    f_dyn.fxx=fxx; f_dyn.fxu=fxu; f_dyn.fuu=fuu;
end

% generate obstacles for obs avoidance problems
    % [h,obs_loc]=obs_generator_2d(num_obs,max_x,max_y,max_rad); 
    [h,obs_loc,obstacles]=obs_manual1;
%     [h,obs_loc,obstacles]=obs_manual0;
%     [h,obs_loc,obstacles]=obs_manual3;
%     [h]=obs_manual2;
%     [h,obs_info]=obs_manual_3D;
%     [h, cart_limit,n_bas]=inv_pend_constraints;
    % [h,w_max, gamma_gs,theta_max, delta_max, Tmin, Tmax]=rocket_const(xf);
%     [h, cart_limit]=cart_pole_constraints;
if barrier_state_button == 1    
% generate the DBaS given h and f
    [f_w,fx_w,fu_w,w0,wf,fxx_w,fxu_w,fuu_w]=DBaS_dyn(x,u,x0,xf,f_dyn,h,ddp_2nd_order);
    % [f_w,fx_w,fu_w,w0,wf]=BaS_dyn(x,u,x0,xf,f,h,n_bas,dt);
    % alpha=0.1;
    % [f_w,fx_w,fu_w,w0,wf]=DBaS_dyn_with_memory(x,u,x0,xf,f,h,alpha,n_bas);
% augment the bas
    [fbar,fbarx,fbaru,fbarxx,fbarxu,fbaruu,xbar0,xbarf,nbar]=Safety_Embedding_dynamics(x,u,x0,xf,f_dyn,w0,wf,f_w,fx_w,fu_w,fxx_w,fxu_w,fuu_w,1,ddp_2nd_order);
% redfine the variables to use vanilla DDP
    f_dyn.f=fbar; f_dyn.fx=fbarx; f_dyn.fu=fbaru;
    if ddp_2nd_order
        f_dyn.fxx=fbarxx; f_dyn.fxu=fbarxu; f_dyn.fuu=fbaruu;
    end
    x0=xbar0; xf=xbarf;n=nbar;
end

sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
%% Quadratic costs (running cost and terminal cost)
% quad :
Q=0*eye(n); 
% Q(n,n)=.2;
% R=1*1e-2*eye(m); 
R=1/2*1e-3*eye(m); 
S=.1*eye(n);
% S(10:12,10:12)=10*eye(3);
% inv pend :
% Q_bf =800;
% Q_bf =0.2;
Q_bf =1e-4;
% R=0.5*1e-2*eye(m);
% S=10*eye(n);
% S(n,n)=0.05;
% S(2,2)=10;
% S=10*eye(n);
% S(1,1)=50;
% S(2,2)=800;
% S(2,2)=500;

S=100*eye(n);
S_bf=0.05;

run_cost=@(x,u,deriv_bool) run_quad_cost(x,u,Q,R,xf,deriv_bool);

% for penalty method comparisons
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
% input constranits
u_min=[];u_max=[];
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
lambda = 1;              % initial value for lambda for regularization
dlambda= 1;              % initial value for dlambda for regularization
lambdaFactor = 1.6;      % lambda scaling factor
lambdaMax = 1e10;        % lambda maximum value
lambdaMin = 1e-6;        % below this value lambda = 0
opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda', dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax, 'lambdaMin', lambdaMin);
%% Discrete DDP
% [X,U,J,lambda,dlambda,alpha,k_u,K_u,ii,iter_succ,L] =disc_ddp_alg(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par);
[X,U,J,lambda,dlambda,alpha,k_u,K_u,ii,iter_succ] =disc_ddp_alg_penalty(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par,h);
%% simulate outside
% plot3(X(10,:),X(11,:),X(12,:),'LineWidth',1); hold on;
% [X,U]=simulate_ddp(f_dyn,k_u,K_u,x0,N,X,U);
% [X_sim,U_sim]=simulate_cbf_ddp(f_dyn,k_u,K_u,x0,N,X,U);
% plottings
% cart_pole_plot(T,X_sim,U_sim,x0,xf,cart_limit);
% cart_pole_plot(T,X,U,x0,xf,cart_limit);
car_2d_plot(T,X,x0,xf,obstacles);
% inv_pend_plot(T,X);
% quad_plot(f_dyn,k_u,K_u,N,U,T,X,x0,xf,obs_info);
% rocket_plot(T,X);

X(:,end)
ii
iter_succ
J(end)
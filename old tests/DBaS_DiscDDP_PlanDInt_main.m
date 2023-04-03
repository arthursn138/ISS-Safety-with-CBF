%% Discrete Time DDP for planar double integrator with obstacles with DISCRETE BaS
clear all; 
close all ; clc
%% add libraries/folders
cost_path = '../costs/';
addpath(genpath(cost_path))                               

ddp_path = '../ddp_algorithm/';
addpath(ddp_path)    
%%
% system dynamics: xk+1=xk+f(xk,uk)*dt
% deriv_boolean if derivative wanted should be true
% [f_dyn]=@(x,u,deriv_bool) SafeEmbd_PlanDInt_dynamics(x,u,deriv_bool,f_dyn,f_z,fx_z,fu_z);
% [f,fx,fu]=f_dyn(x,u,deriv_bool);
n=5; m=2;               % states and input dimensions resp
% sampling rate and no. samples
dt = 0.01;
N = 800;
% initial states
x0=zeros(4,1);
% x0=[1;5;0;0];
% terminal states (desired state)
xf=zeros(4,1);
xf=[4;4;0;0];
% define safety h= (x1-o1)^2+(x2-o2)^2-r^2 
o11=2.5; o12=2; r1=1/2; circ1 = [o11-r1 o12-r1 2*r1 2*r1];
o21=2; o22=4; r2=1/2; circ2 = [o21-r2 o22-r2 2*r2 2*r2];
o31=1; o32=2.5; r3=1/2; circ3 = [o31-r3 o32-r3 2*r3 2*r3];
o41=4.5; o42=1; r4=1/2; circ4 = [o41-r4 o42-r4 2*r4 2*r4];
safe_par = struct('o11',o11,'o12',o12,'r1',r1,'o21',o21,'o22',o22,'r2',r2,'o31',o31,'o32',o32,'r3',r3,'o41',o41,'o42',o42,'r4',r4);
% planar dyn
[f]=@(x,u) PlanDInt_dynamics(x,u,dt,0);
% symbolic DBaS
[f_z,fx_z,fu_z,z0,zf]=DBaS_dyn_symb(x0,xf,f,safe_par);

% new ic and fc
x0=[x0;z0];
xf=[xf;zf];

[f_dyn]=@(x,u,deriv_bool) SafeEmbd_PlanDInt_dynamics(x,u,deriv_bool,f_z,fx_z,fu_z,dt);

ddp_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
%% Quadratic costs (running cost and terminal cost)
Q=1e-3*eye(n); Q(5,5)=1e-3; R=.01*eye(m); S=100*eye(n);% S(5,5)=0.1; %state, input and term.cond matrices
run_cost=@(x,u,deriv_bool) run_quad_cost(x,u,Q,R,xf,deriv_bool);
% [L,L_x, L_u, L_xu, L_uu, L_xx]=run_cost(x,u,deriv_bool);
term_cost=@(x,deriv_bool) terminal_quad_cost(x,xf,S,deriv_bool);
% [l_term,l_term_x,l_term_xx] = term_cost(x,xf,S,deriv_bool)
% input constranits
u_min=[];u_max=[];
% nominal input and state
ubar = 0.1*ones(m, N-1);   % nominal control
xbar(:,1) = x0;  % initial state
for k=1:N-1
    xbar(:,k + 1) = f_dyn(xbar(:, k), ubar(:, k),false);
end 
%
iter= 20;           % number of max iterations
toler= 1e-3;              % cost change 1e-3
% for regularization:
lambda = 1;          % initial value for lambda for regularization
dlambda= 1;         % initial value for dlambda for regularization
lambdaFactor = 1.6;      % lambda scaling factor
lambdaMax = 1e10;        % lambda maximum value
lambdaMin = 1e-6;        % below this value lambda = 0
opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda', dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax, 'lambdaMin', lambdaMin);
% Discrete DDP function
[X,U,J,lambda,dlambda,alpha] =disc_ddp_alg(f_dyn,run_cost,term_cost,ddp_par,ubar,xbar,opt_par);

% plotting
T=0:dt:dt*N-dt;
figure(1)
subplot(2,2,1);
plot(xf(1),xf(2),'x'); 
hold on
scatter(X(1,1),X(2,1),'o'); 
rectangle('Position',circ1,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
rectangle('Position',circ2,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
rectangle('Position',circ3,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
rectangle('Position',circ4,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
plot(X(1,:),X(2,:)); 
title('x_1 vs x_2')
subplot(2,2,2)
plot(T,X(1,:),'--',T,X(2,:),':'); title('x_1, x_2 vs t'); legend('x_1','x_2');
subplot(2,2,3)
plot(T(2:end),U(1,:),'--',T(2:end),U(2,:),':'); title('u_1, u_2 vs t'); legend('u_1','u_2');
subplot(2,2,4)
plot(J); title('Cost vs iteration');
figure(2)
plot(T,X(5,:));
title('BaS vs time')

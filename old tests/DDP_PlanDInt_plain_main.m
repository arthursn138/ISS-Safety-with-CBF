%% Discrete Time DDP for planar double integrator:
clear all; close all ; clc

% system dynamics: xk+1=xk+f(xk,uk)*dt
% deriv_boolean if derivative wanted should be true
[f_dyn]=@(x,u,deriv_bool) PlanDInt_dynamics(x,u,deriv_bool);
% [f,fx,fu]=f_dyn(x,u,deriv_bool);
n=4; m=2;               % states and input dimensions resp
% sampling time and no. samples
dt = 0.02;
N = 400;
% initial states
x0=zeros(n,1);
x0=[4;4;0;0];
% terminal states (desired state)
xf=[0;0;0;0];

ddp_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
% Quadratic costs (running cost and terminal cost)
Q=eye(4); R=eye(2); S=10*eye(4); %state, input and term.cond matrices
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
    xbar(:,k + 1) = xbar(:, k)+ dt*f_dyn(xbar(:, k), ubar(:, k),false);
end 
% optimization parameters
% unconstrained case cooment rho_x and rho_u
% Opt.rho_x = 1;
% Opt.rho_u = 1;
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

%% plotting
T=0:dt:dt*N-dt;
figure(2)
subplot(2,2,1);plot(xf(1),xf(2),'x'); 
hold on
plot(X(1,1),X(2,1),'o'); hold on;
plot(X(1,:),X(2,:)); 
title('x_1 vs x_2')
subplot(2,2,2)
plot(T,X(1,:),'--',T,X(2,:),':'); title('x_1, x_2 vs t'); legend('x_1','x_2');
subplot(2,2,3)
plot(T(2:end),U(1,:),'--',T(2:end),U(2,:),':'); title('u_1, u_2 vs t'); legend('u_1','u_2');
subplot(2,2,4)
plot(J); title('Cost vs iteration');



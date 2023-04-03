%% Discrete Time DDP for planar double integrator with obstacles wiht CONTINUOUS BaS 
clear all; 
% close all ; clc

% system dynamics: xk+1=xk+f(xk,uk)*dt
% deriv_boolean if derivative wanted should be true
% [f_dyn]=@(x,u,deriv_bool) SafeEmbd_PlanDInt_dynamics(x,u,deriv_bool,f_dyn,f_z,fx_z,fu_z);
% [f,fx,fu]=f_dyn(x,u,deriv_bool);
n=5; m=2;               % states and input dimensions resp
% sampling rate and no. samples
dt = 0.01;
N = 400;
% initial states
x0=zeros(4,1);
x0=[1;5;0;0];
% terminal states (desired state)
xf=zeros(4,1);
% xf=[3;3;0;0];
% define safety h= (x1-o1)^2+(x2-o2)^2-r^2 
o11=0; o12=1.75; r1=1/2; pos1 = [o11-r1 o12-r1 2*r1 2*r1];
o21=.75; o22=1.75; r2=1/2; pos2 = [o21-r2 o22-r2 2*r2 2*r2];
o31=1.5; o32=1.6; r3=1/2; pos3 = [o31-r3 o32-r3 2*r3 2*r3];
o41=1; o42=-0.5; r4=1/2; pos4 = [o41-r4 o42-r4 2*r4 2*r4];
gamma=0;  
% symbolic BaS
[f_z,fx_z,fu_z,z0,zf]=BaS_dyn_symb(o11,o12,r1,o21,o22,r2,o31,o32,r3,o41,o42,r4,gamma,x0,xf);
% symbolic to function handle
f_z=matlabFunction(f_z);
fx_z=matlabFunction(fx_z);
fu_z=matlabFunction(fu_z);

% new ic and fc
x0=[x0;z0];
xf=[xf;zf];

[f_dyn]=@(x,u,deriv_bool) SafeEmbd_PlanDInt_dynamics(x,u,deriv_bool,f_z,fx_z,fu_z);

ddp_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
%% Quadratic costs (running cost and terminal cost)
Q=1e-3*eye(n); Q(5,5)=1e-6; R=.01*eye(m); S=100*eye(n);% S(5,5)=0.1; %state, input and term.cond matrices
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
rectangle('Position',pos1,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
rectangle('Position',pos2,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
rectangle('Position',pos3,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
rectangle('Position',pos4,'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
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

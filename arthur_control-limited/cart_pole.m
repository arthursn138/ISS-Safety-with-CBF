%% Control-Limited Differential Dynamic Programming
% Arthur Nascimento, Hassan Almubarak
% nascimento@gatech.edu, halmubarak@gatech.edu
% ACDS Lab @ Georgia Tech
% Last Update March/2022

% This is an example for Cart-Pole system with Barrier States and control-
% limited discrete DDP.
% Follow those steps to run any example
% 1. Call the system's dynamics
% 2. Define control limits: first collumn is the lower bound, second
% collumn is the upper bound. One row per control (u_lim should be m-sized)
% To skip this step, simply leave u_lims = []
% 3. Call DBaS_dyn to generate the DBaS dynamics
% 4. Call Safety_Embedding_dynamics to augment the DBaS to the system's
% dynamics
% 5. Define numerical (DDP and optimization) paramters and run the discrete
% control-limited DDP.

% To be implemented: graphically generate obstacles

clear; close all; clc;

% Add Paths
parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory)); 

%% Initialize system parameters

% Numerical parameters
dt = 0.02;  % Sampling
N = 250;    % Horizon (# of iteration that discretizes each pass)
% T = dt*N;
T = 0:dt:dt*N-1*dt;
% T = [0:dt:dt*N*dt dt*NT = 0:dt:dt*N*dt];

% Set initial and final state
x0=[0;0;0;0];
xf=[0;pi;0;0];

% System Dynamics
[f,fx,fu,fxx,fxu,fuu,x,u] = cart_pole_dynamics(dt,1);

n = length(x); % State dimensions
m = length(u); % Input dimensions

if length(xf)~=n || length(x0)~=n
    error('wrong dimention of boundary conditions');
end

ddp_2nd_order = 0; % 1 For Full DDP (not fully implemented on this code)

f_dyn.f = f; f_dyn.fx = fx; f_dyn.fu = fu; 
f_dyn.fxx = fxx; f_dyn.fxu = fxu; f_dyn.fuu = fuu;

% ====== Barrier states (set to 0 to skip BaS) =======
barrier_state_button = 1;

% Generate obstacles and safety condition
[h, cart_limit] = cart_pole_constraints;

if barrier_state_button == 1
% Generate the Discrete Barrier States given h and f
    [f_w,fx_w,fu_w,w0,wf,fxx_w,fxu_w,fuu_w] = DBaS_dyn(x,u,x0,xf,f_dyn,...
        h,ddp_2nd_order);
% Augment the BaS 
    [fbar,fbarx,fbaru,fbarxx,fbarxu,fbaruu,xbar0,xbarf,nbar]=Safety_Embedding_dynamics(x,...
        u,x0,xf,f_dyn,w0,wf,f_w,fx_w,fu_w,fxx_w,fxu_w,fuu_w,1,ddp_2nd_order);
% Redfine the variables to use DDP
    f_dyn.f = fbar; f_dyn.fx = fbarx; f_dyn.fu = fbaru;
    x0 = xbar0; xf = xbarf; n = nbar;
end

sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf, 'n', n, 'm', m);

%% Quadratic costs (running cost and terminal cost)

Q = 0*eye(n);   % States weights
Q(1,1) = 0.1;   % States weights
Q(2,2) = 0.2;   % States weights
% Q(3,3) = 1;   % States weights
% Q(4,4) = 1;   % States weights
R = 0.2*1e-2*eye(m);    % Controls weights
S = 1000*eye(n);        % Final cost
S(3,3) = 100;           % Final cost
S(4,4) = 500;           % Final cost

% % % quad :
% % Q=0*eye(n); 
% % R=0.5*1e-1*eye(m); 
% % S=10*eye(n);
% % S(1,1)=50;
% % S(2,2)=800;
% % % S(2,2)=500;

Q_bf = 1;
S_bf = 0.05;

% Cost functions
run_cost = @(x,u,deriv_bool) run_quad_cost(x, u, Q, R, xf, deriv_bool);
if barrier_state_button == 1
    Q(n,n)=Q_bf;
    S(n,n)=S_bf;
    run_cost=@(x,u,deriv_bool) run_quad_cost(x, u, Q, R, xf, deriv_bool);
end
term_cost = @(x,deriv_bool) terminal_quad_cost(x, xf, S, deriv_bool);

% Input constraints
u_lims  = [-10 20];  % Force on the x direction (N), positive to the right
% u_lims = [];          % Use this if you don't want to limit your inputs

[r,c] = size(u_lims);
if ~isempty(u_lims) && r ~= m
    error('Wrong dimension of control limits (# of rows)');
elseif ~isempty(u_lims) && c ~= 2
    error('Wrong dimension of control limits. # of collumns should be 2: first for lower bound, second for upper bound');
end

%% Nominal input and state

ubar = 0.0*ones(m, N-1);    % nominal control
xbar = []; 
xbar(:,1) = x0;             % initial state

% Takes initial condition and nominal control and propagates the states 
% according to the system's dynamics
for k=1:N-1
    xbar(:,k + 1) = f_dyn.f(xbar(:,k), ubar(:,k));
end

%% Optimization parameters

iter = 500;             % Maximum number of iterations
toler = 1e-3;           % For cost change

% For the regularization part by Yuichiro Aoyama
lambda = 0.5;            % Initial value for lambda for regularization
dlambda= 1;              % Initial value for dlambda for regularization
lambdaFactor = 1.6;      % Lambda scaling factor
lambdaMax = 1e10;        % Lambda maximum value
lambdaMin = 1e-6;        % Below this value lambda = 0

opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda',...
    dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax,...
    'lambdaMin', lambdaMin);

%% Discrete Control-Limited DDP

[X,U,J,lambda,dlambda,alpha,k_u,K_u,ii,iter_succ,L] = ctrl_lim_disc_ddp_alg(ddp_2nd_order,...
    f_dyn, run_cost, term_cost, sys_par, ubar, xbar, u_lims, opt_par, h);

%% Checking final values
X(:,end);
ii;
iter_succ;
J(end)

%% Plottings
cart_pole_plot(T,X,U,x0,xf,u_lims,cart_limit);


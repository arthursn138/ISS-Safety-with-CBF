%% Control-Limited Differential Dynamic Programming
% Arthur Nascimento, Hassan Almubarak
% nascimento@gatech.edu, halmubarak@gatech.edu
% ACDS Lab @ Georgia Tech
% Last Update March/2022

% This is an example for 2D single integrator.
% Follow those steps to run any example
% 1. call the system's dynamics
% 2. Define DDP and optimization paramters and run vanilla DDP

% To be implemented: generate obstacles, define safe set functionand call 
% DBaS_dyn and Safety_Embedding_dynamics

clear; close all; clc;

% Add Paths
parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory)); 

%% Initialize system parameters

% Numerical parameters
dt = 0.02;  % Sampling
N = 150;    % Horizon (# of iteration that discretizes each pass)
T = 0:dt:dt*N-1*dt;
% T = dt*N;

% Set initial and final state
x0 = [-1;3];
xf = [1;3];

% System Dynamics
[f,fx,fu,fxx,fxu,fuu,x,u] = single_integrator_dynamics(dt,1);

n = length(x); % State dimensions
m = length(u); % Input dimensions

f_dyn.f = f; f_dyn.fx = fx; f_dyn.fu = fu; 
f_dyn.fxx = fxx; f_dyn.fxu = fxu; f_dyn.fuu = fuu;

sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf, 'n', n, 'm', m);

%% Quadratic costs (running cost and terminal cost)

Q = 0*eye(n);       % States weights
R = 1e-2*eye(m);    % Controls weights
S = 10*eye(n);      % Final cost

% Cost functions
run_cost = @(x,u,deriv_bool) run_quad_cost(x, u, Q, R, xf, deriv_bool);
term_cost = @(x,deriv_bool) terminal_quad_cost(x, xf, S, deriv_bool);

% Input constraints
u_lims  = [-.5 .5;      % controller 1 limits (no units for this model)
             -2  2];    % controller 2 limits (no units for this model)
% u_lims = [];

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

%% Discrete DDP

ddp_2nd_order = 0; % For Full DDP
[X,U,J,lambda,dlambda,alpha,k_u,K_u,ii,iter_succ,L] = ctrl_lim_disc_ddp_alg(ddp_2nd_order,...
    f_dyn, run_cost, term_cost, sys_par, ubar, xbar, u_lims, opt_par);

%% Checking final values
X(:,end)
% ii
% iter_succ
J(end)

%% Plottings
car_2d_plot(T,X,U,x0,xf);


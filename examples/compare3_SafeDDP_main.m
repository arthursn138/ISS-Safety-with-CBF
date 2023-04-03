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
% clear
% close all
% clc

%% add paths
% cost_path = '../costs/';
% addpath(genpath(cost_path))                               

% ddp_path = '../ddp_algorithm/';
% addpath(ddp_path)     

% sys_path = '../Systems Dynamics/';
% addpath(sys_path)     

% genral_path = '../hassan_safety_embedded_ddp/';
% addpath(genpath(genral_path))  
parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory)); 
%% initialize system parameters
dt=0.02; %sampling
N=300; % horizon
T=dt*N;
T=0:dt:dt*N-1*dt;

[f,fx,fu,fxx,fxu,fuu,x,u]=Diff_wheel_dynamics(dt,1);

f_dyn.f=f; f_dyn.fx=fx; f_dyn.fu=fu; 
f_dyn.fxx=fxx; f_dyn.fxu=fxu; f_dyn.fuu=fuu;
ddp_2nd_order = 1;


counter_1_penalty=0;
counter_1ov2_penalty=0;
counter_1ov4_penalty=0;
counter_1ov10_penalty=0;
counter_1_bas=0;
counter_1ov2_bas=0;
counter_1ov4_bas=0;
counter_1ov10_bas=0;
counter_flags=0;
XX1=[];XX2=[];XX3=[];XX4=[];

num_trials=10; 

NUM_OBS=[1 2 3 4 5 6 7 8 9 10];
NUM_OBS=3;

BsS_exp_1ov4=[];
Penalty_exp_1ov4=[];

BsS_exp_1ov10=[];
Penalty_exp_1ov10=[];
JJ_BaS_means=[];
JJ_Penalty_means=[];
X1_bas=[];
X1_pen=[];
X2_bas=[];
X2_pen=[];
ALLOBS=[];
X0=[];Xf=[];
ITER_PEN=[];
ITER_BAS=[];
ITER_SUCC_BAS=[];
ITER_SUCC_PEN=[];

CPUTIME_PEN=[];
CPUTIME_BAS=[];

for nn=1:length(NUM_OBS);
    nn
    num_obs=NUM_OBS(nn);
    
    counter_1_penalty=0;
    counter_1ov2_penalty=0;
    counter_1ov4_penalty=0;
    counter_1ov10_penalty=0;
    counter_1_bas=0;
    counter_1ov2_bas=0;
    counter_1ov4_bas=0;
    counter_1ov10_bas=0;
    counter_flags=0;

    JJ_BaS_nobs=[];
    JJ_Penalty_nobs=[];
    iter_bas=[];
    iter_succ_bas=[];
    iter_pen=[];
    iter_succ_pen=[];
    cputime_pen=[];
    cputime_bas=[];
    
for ll=1:num_trials
%     ll
    
    barrier_state_button = 0;
    penalty_button = 1;

    [h,obs_loc,obstacles]=obs_random(num_obs);
    n=length(x); %state dim (total, with Bas)
    m=length(u); %input dim
    % specify initial state for the original system
    x0 = [-3;0;0]+(rand(3,1)-0.5);
    X0=[X0 x0];
    
    xf=[3;0;0]+(rand(3,1)-0.5);
    Xf=[Xf xf];

    flag_obs_cover_ic=0;
%         for kk=1:length(h)
%             if ( (h{kk}(x0) <0) || (h{kk}(xf)) <0)
%             flag_obs_cover_ic=1;
%             break;
%             end
%         end
    kk=0;
    while kk<length(h)
          kk=kk+1;
          if ( (h{kk}(x0) <0) || (h{kk}(xf)) <0)
             [h,obs_loc,obstacles]=obs_random(num_obs);
             k=0;
          end
    end
    ALLOBS=[ALLOBS;obstacles];
    
        if flag_obs_cover_ic==1
            counter_flags=counter_flags+1;
            continue;
        end
    
    f_dyn.f=f; f_dyn.fx=fx; f_dyn.fu=fu; 
    f_dyn.fxx=fxx; f_dyn.fxu=fxu; f_dyn.fuu=fuu;
    sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf,'n',n,'m',m);
%% Quadratic costs (running cost and terminal cost)
% quad :
    Q=0*eye(n); 
% Q(n,n)=.2;
% R=1*1e-4*eye(m); 
% S=1*eye(n);
% S(10:12,10:12)=2*eye(3);
% inv pend :
% Q_bf =800;
    Q_bf =1e-5;
    R=0.5*1e-4*eye(m);
    S=20*eye(n);

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

    term_cost=@(x,deriv_bool) terminal_quad_cost(x,xf,S,deriv_bool);
    % input constranits
    u_min=[];u_max=[];
%% nominal input and state
    ubar = 0*ones(m, N-1); % nominal control
    xbar=[]; xbar(:,1) = x0;          % initial state
    for k=1:N-1
        xbar(:,k + 1) = f_dyn.f(xbar(:, k), ubar(:, k));
    end
%% optimization parameters
    iter= 200;                % number of max iterations
    toler= 1e-3;             % cost change 1e-3
% for regularization part by Yuichiro Aoyama
    lambda = 1;              % initial value for lambda for regularization
    dlambda= 1;              % initial value for dlambda for regularization
    lambdaFactor = 1.6;      % lambda scaling factor
    lambdaMax = 1e10;        % lambda maximum value
    lambdaMin = 1e-6;        % below this value lambda = 0
    opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda', dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax, 'lambdaMin', lambdaMin);
%% Discrete DDP
% [X,U,J,lambda,dlambda,alpha,k_u,K_u] =disc_ddp_alg(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par);

tic
[X,U,J,lambda,dlambda,alpha,k_u,K_u,iter,iter_succ] =disc_ddp_alg_penalty(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par,h);
cur_cputime=toc;
cputime_pen=[cputime_pen cur_cputime];
X1_pen=[X1_pen;X(1,:)];
X2_pen=[X2_pen;X(2,:)];
iter_pen=[iter_pen iter];
if isempty(iter_succ)
    iter_succ=0;
end

iter_succ_pen=[iter_succ_pen iter_succ];

%     car_2d_plot(T,X,x0,xf,obstacles);

    L=0; %initialze
    for jj=1:N-1
        L=L+run_cost(X(:,jj),U(:,jj),false);
    end
    [L_f]= term_cost(X(:,end),false);
    cost=L+L_f;
    JJ_Penalty_nobs=[JJ_Penalty_nobs cost];

    %% 
    xf_real=X(:,end);
    if norm(xf_real(1:2)-xf(1:2)) <= 1
        counter_1_penalty=counter_1_penalty+1;
    end
    if norm(xf_real(1:2)-xf(1:2)) <= 0.5
        counter_1ov2_penalty=counter_1ov2_penalty+1;
    end
    if norm(xf_real(1:2)-xf(1:2)) <= 0.25
        counter_1ov4_penalty=counter_1ov4_penalty+1;
    end
    if norm(xf_real(1:2)-xf(1:2)) <= 0.1
        counter_1ov10_penalty=counter_1ov10_penalty+1;
    end
%% bas


    barrier_state_button = 1;
    penalty_button = 0;

% generate obstacles for obs avoidance problems
    if barrier_state_button == 1    
% generate the DBaS given h and f
    [f_w,fx_w,fu_w,w0,wf,fxx_w,fxu_w,fuu_w]=DBaS_dyn(x,u,x0,xf,f_dyn,h,ddp_2nd_order);
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
% R=1*1e-4*eye(m); 
% S=1*eye(n);
% S(10:12,10:12)=2*eye(3);
% inv pend :
% Q_bf =800;
Q_bf =1e-5;
R=0.5*1e-4*eye(m);
S=20*eye(n);
if barrier_state_button == 1    
S(n,n)=0.05;
end

if barrier_state_button == 1    
    Q(n,n)=Q_bf;
    run_cost=@(x,u,deriv_bool) run_quad_cost(x,u,Q,R,xf,deriv_bool);
end

term_cost=@(x,deriv_bool) terminal_quad_cost(x,xf,S,deriv_bool);
% input constranits
u_min=[];u_max=[];
%% nominal input and state
ubar = 0.0*ones(m, N-1); % nominal control
xbar=[]; xbar(:,1) = x0;          % initial state
for k=1:N-1
    xbar(:,k + 1) = f_dyn.f(xbar(:, k), ubar(:, k));
end
%% optimization parameters
iter= 200;                % number of max iterations
toler= 1e-3;             % cost change 1e-3
% for regularization part by Yuichiro Aoyama
lambda = 1;              % initial value for lambda for regularization
dlambda= 1;              % initial value for dlambda for regularization
lambdaFactor = 1.6;      % lambda scaling factor
lambdaMax = 1e10;        % lambda maximum value
lambdaMin = 1e-6;        % below this value lambda = 0
opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda', dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax, 'lambdaMin', lambdaMin);
%% Discrete DDP
% [X,U,J,lambda,dlambda,alpha,k_u,K_u] =disc_ddp_alg(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par);
tic
[X,U,J,lambda,dlambda,alpha,k_u,K_u,iter,iter_succ] =disc_ddp_alg_penalty(ddp_2nd_order,f_dyn,run_cost,term_cost,sys_par,ubar,xbar,opt_par,h);
cur_cputime=toc;
cputime_bas=[cputime_bas cur_cputime];
% car_2d_plot(T,X,x0,xf,obstacles);
X1_bas=[X1_bas;X(1,:)];
X2_bas=[X2_bas;X(2,:)];
iter_bas=[iter_bas iter];
if isempty(iter_succ)
    iter_succ=0;
end
iter_succ_bas=[iter_succ_bas iter_succ];
% close all
    L=0; %initialze
    for jj=1:N-1
        L=L+run_cost(X(:,jj),U(:,jj),false);
    end
    [L_f]= term_cost(X(:,end),false);
    cost=L+L_f;
    
    JJ_BaS_nobs=[JJ_BaS_nobs cost];

%% 

    xf_real=X(:,end);
    if norm(xf_real(1:2)-xf(1:2)) <= 1
        counter_1_bas=counter_1_bas+1;
    end
    if norm(xf_real(1:2)-xf(1:2)) <= 0.5
        counter_1ov2_bas=counter_1ov2_bas+1;
    end
    if norm(xf_real(1:2)-xf(1:2)) <= 0.25
        counter_1ov4_bas=counter_1ov4_bas+1;
    end
    if norm(xf_real(1:2)-xf(1:2)) <= 0.1
        counter_1ov10_bas=counter_1ov10_bas+1;
    end
    
end
peformed_trials= num_trials - counter_flags

counter_1_penalty;
counter_1ov2_penalty;
counter_1ov4_penalty;
counter_1ov10_penalty;
succ_rate_1_penalty= counter_1_penalty/peformed_trials*100;
succ_rate_1ov2_penalty= counter_1ov2_penalty/peformed_trials*100;
succ_rate_1ov4_penalty= counter_1ov4_penalty/peformed_trials*100;
succ_rate_1ov10_penalty= counter_1ov10_penalty/peformed_trials*100;


counter_1_bas;
counter_1ov2_bas;
counter_1ov4_bas;
counter_1ov10_bas;
succ_rate_1_bas= counter_1_bas/peformed_trials*100;
succ_rate_1ov2_bas= counter_1ov2_bas/peformed_trials*100;
succ_rate_1ov4_bas= counter_1ov4_bas/peformed_trials*100;
succ_rate_1ov10_bas= counter_1ov10_bas/peformed_trials*100;

BsS_exp_1ov4=[BsS_exp_1ov4 succ_rate_1ov4_bas];
BsS_exp_1ov10=[BsS_exp_1ov10 succ_rate_1ov10_bas];
Penalty_exp_1ov4=[Penalty_exp_1ov4 succ_rate_1ov4_penalty];
Penalty_exp_1ov10=[Penalty_exp_1ov10 succ_rate_1ov10_penalty];

JJ_BaS_nobs_mean=sum(JJ_BaS_nobs)/peformed_trials;
JJ_BaS_means=[JJ_BaS_means JJ_BaS_nobs_mean];
JJ_Penalty_nobs_mean=sum(JJ_Penalty_nobs)/peformed_trials;
JJ_Penalty_means=[JJ_Penalty_means JJ_Penalty_nobs_mean];

ITER_PEN=[ITER_PEN;iter_pen];
ITER_SUCC_PEN=[ITER_SUCC_PEN;iter_succ_pen];
ITER_BAS=[ITER_BAS;iter_bas];
ITER_SUCC_BAS=[ITER_SUCC_BAS;iter_succ_bas];

CPUTIME_BAS=[CPUTIME_BAS;cputime_bas];
CPUTIME_PEN=[CPUTIME_PEN;cputime_pen];

end

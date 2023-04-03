function [x, u,J,k_u,K_u,ii,iter_succ,total_cost] =simple_ddp_alg(f_dyn,run_cost,term_cost,ddp_par,ubar,xbar,opt_par);
% first order discrete ddp (iLQR
    dt=ddp_par.dt; N=ddp_par.N; n=ddp_par.n; m=ddp_par.m; 
    x0=ddp_par.x0; xf=ddp_par.xf;
    iter=opt_par.iter; toler=opt_par.toler;
    conv_flag=0;
    % compute nominal running cost
    L=0; %initialze
    for jj=1:N-1
        L=L+run_cost(xbar(:,jj),ubar(:,jj),false);
    end
    % compute nominal terminal cost
    [L_f]= term_cost(xbar(:,end),false);
    total_cost=L+L_f;
    % initialize the costs for eah iteration
    J=zeros(1, iter);
    %limits of u
    lims=[];
    % initialize expected decrease in cost
    dJ_exp=[0 0];
    iter_succ=[];
    % discrete DDP algorithm
    for ii=1:iter
        J(ii)=total_cost;                
        % backward pass: from N to 1
        % for N:
        % get V,Vx,Vxx
        [L_f,L_f_x,L_f_xx]= term_cost(xbar(:,end),true);
        V_plus=L_f; Vx_plus=L_f_x; Vxx_plus=L_f_xx;
        % initialize ff gain and fb gain
        k_u = zeros(m,N); % ff gain
        K_u = zeros(m,n,N); % fb gain
        for jj=N-1:-1:1
            f_xk=f_dyn.fx(xbar(:,jj),ubar(:,jj));
            f_uk=f_dyn.fu(xbar(:,jj),ubar(:,jj));
            [~,L_xk,L_uk,L_xuk,L_xxk,L_uuk] = run_cost(xbar(:,jj), ubar(:,jj), true);
            Qx=L_xk + f_xk'*Vx_plus;
            Qu=L_uk + f_uk'*Vx_plus;
            
            Qxx=L_xxk + f_xk'*Vxx_plus*f_xk;
            Qxu=L_xuk + f_xk'*Vxx_plus*f_uk; Qux = Qxu';
            Quu=L_uuk + f_uk'*Vxx_plus*f_uk;
            
            
            % making sure Quu is not NAN then not negative def
            if any(isnan(Quu))
                error('Quu_reg has NaN in ddp\n')
            end
            if min(eig(Quu))<0
                fprintf(1, 'Quu ndf at %d . make regularizer large \n',jj);
                break            
            end
            
            
            k_u(:,jj)= -Quu\Qu; % ff gain
            K_u(:,:,jj)= -Quu\Qux;  % fb gain
            
            Vx_plus = Qx + K_u(:,:,jj)'*Quu*k_u(:,jj) + K_u(:,:,jj)'*Qu + Qxu*k_u(:,jj); %Qx + Qxu*l(:,k);
            Vxx_plus = Qxx + K_u(:,:,jj)'*Quu*K_u(:,:,jj) + K_u(:,:,jj)'*Qux + Qxu*K_u(:,:,jj);
            
            % making sure Vxx is positive definite (numerical issues)
            Vxx_plus = 0.5 * (Vxx_plus + Vxx_plus'); 
            % for expected decrease in cost:
            dJ_exp = dJ_exp + [k_u(:,jj)'*Qu  0.5*k_u(:,jj)'*Quu*k_u(:,jj)];
        end
        
        % forward pass:
        % initilizations 
        x = zeros(n,N);
        x(:,1) = x0;
        u = zeros(m,N-1);
        deltax = x(:,1) - xbar(:,1);
        % start forward propagation
        L = 0;
        for jj=1:N -1
            deltau = k_u(:,jj) + K_u(:,:,jj)*deltax;
            u(:,jj) = ubar(:,jj) + deltau;
            x(:,jj + 1)=  f_dyn.f(x(:,jj),u(:,jj));          
            deltax = x(:,jj + 1) - xbar(:,jj + 1);  
            L = L + run_cost(x(:,jj),u(:,jj),false);
        end
        L_f = term_cost(x(:,end),true);
        total_cost_new = L + L_f;
        % actual decrease in cost
        dcost = total_cost-total_cost_new;
               
        if abs(dcost) < toler           
           conv_flag = 1;
        end
        
        % overwriting for the next iteration
        total_cost=total_cost_new;          
        xbar = x; ubar = u;
        if conv_flag, fprintf('convergence of DDP dcost: %.3f\n',dcost);
           break; 
        end
        if isinf(total_cost)
                fprintf('diverged, try with larger lambda: %2f\n',lambda)
        end
    end
    iter_succ=min(iter_succ);J = J(1:ii);
end
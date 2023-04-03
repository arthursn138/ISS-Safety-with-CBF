%% get the control gains (ff and fb) and simulate the system:
function [x,u]=simulate_quadrotor(f_dyn,k_u,K_u,x0,N,X,U,dt)
    x(:,1)=x0; xbar= X; ubar =U;
    dt=0.01;
    t=0:dt:dt*N-1*dt;
    sinoft = sin(t); amp=15;
    for jj=1:N-1
%         if (x(1,jj)>1.5)
%             x(1,jj)=1.5;
%         else if (x(1,jj)<-1.5)
%             x(1,jj)=-1.5;
%             end
%         end
        
        deltax = x(:,jj) - xbar(:,jj);   
        deltau = 0 + K_u(:,:,jj)*deltax;        
        u(:,jj) = ubar(:,jj) + deltau;
        
        x(:,jj + 1)= f_dyn.f(x(:,jj),u(:,jj));

        x(7,jj+1) = x(7,jj+1) + dt*amp*randn*sinoft(jj);
        x(8,jj+1) = x(8,jj+1) + dt*amp*randn*sinoft(jj);
        x(9,jj+1) = x(9,jj+1) + dt*amp*randn*sinoft(jj);

    end
end
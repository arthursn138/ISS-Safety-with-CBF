%% get the control gains (ff and fb) and simulate the system:
function [x,u]=simulate_cbf_ddp(f_dyn,k_u,K_u,x0,N,X,U)
    x(:,1)=x0; xbar= X; ubar =U;
    options = optimoptions('fmincon','Display','off');
    % systems paramters
    mp = 0.05;      % pole mass (kg)
    mc = 1;         % cart mass (kg)
    g = 9.8;        % gravitational force (m/s^2)
    l = 2;          % pole length (m)
%     alpha1=1000;
%     alpha2=200;
%     alpha1=500;
%     alpha2=100;
%     alpha1=50;
%     alpha2=10;
    alpha1=250;
    alpha2=25;
% best response:
    alpha1=100;
    alpha2=20;
    xlim=1.5;
    w1=0;w2=0;W1=[];W2=[];
    for jj=1:N-1
        deltax = x(:,jj) - xbar(:,jj);   
        deltau = K_u(:,:,jj)*deltax;        
        u(:,jj) = ubar(:,jj) + deltau;
        %filter with CBF:    
        unom=u(:,jj);
%         [u_cbf]=fmincon(@(u_cbf) cbf_objective(u_cbf,unom),[unom;w1;w2],[],[],[],[],[-30;-10000;-10000],[30;10000;10000],@(u_cbf) cbf_constraint(u_cbf,x(:,jj),mp,l,g,mc,alpha,xlim),options);
%         [u_cbf]=fmincon(@(u_cbf) cbf_objective(u_cbf,unom),[unom;w1;w2],[],[],[],[],[],[],@(u_cbf) cbf_constraint(u_cbf,x(:,jj),mp,l,g,mc,xlim),options);
%         [u_cbf]=fmincon(@(u_cbf) cbf_objective(u_cbf,unom),[unom],[],[],[],[],[-15],[15],@(u_cbf) cbf_constraint(u_cbf,x(:,jj),mp,l,g,mc,alpha1,alpha2,xlim),options);
        [u_cbf]=fmincon(@(u_cbf) cbf_objective(u_cbf,unom),[unom],[],[],[],[],[],[],@(u_cbf) cbf_constraint(u_cbf,x(:,jj),mp,l,g,mc,alpha1,alpha2,xlim),options);
        u(:,jj) =u_cbf(1);
%         w1=u_cbf(2);
%         w2=u_cbf(3);
        W1=[W1 w1];
        W2=[W2 w2];
        %apply to system:
        x(:,jj + 1)= f_dyn.f(x(:,jj),u(:,jj));    
    end
    figure(10000)
    plot(W1);hold on; plot(W2); hold off;
    
end
% manual
function [c,ceq] = cbf_constraint(u,x,mp,l,g,mc,alpha1,alpha2,xlim)
c =-(-2*x(1)*((mp*sin(x(2)))*(l*x(4)^2+g*cos(x(2)))+u(1))/(mc+mp*(sin(x(2)))^2)-2*x(3)^2+alpha1*(xlim^2-x(1)^2)+alpha2*-2*x(1)*x(3));
ceq = [];
end
 % CBF Obj function
function f= cbf_objective(u,unom)
f=(u(1)-unom)'*(u(1)-unom);
end

% optimal decay:
 % CBF Obj function
% function f= cbf_objective(u,unom)
% f=(u(1)-unom)'*(u(1)-unom)+20*u(2)^2+20*u(3)^2;
% end
% function [c,ceq] = cbf_constraint(u,x,mp,l,g,mc,xlim)
% c =-(-2*x(1)*((mp*sin(x(2)))*(l*x(4)^2+g*cos(x(2)))+u(1))/(mc+mp*(sin(x(2)))^2)-2*x(3)^2+u(2)*(xlim^2-x(1)^2)+u(3)*-2*x(1)*x(3));
% ceq = [];
% end



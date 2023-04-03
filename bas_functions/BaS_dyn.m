%% Dynamics of Continous Time Barrier States (BaS)
function [f_w,fx_w,fu_w,w0,wf]=BaS_dyn(x,u,x0,xf,f_dyn,h,n_bas,dt)
    xplus=f_dyn(x,u); dx=(xplus-x); xdot=dx/dt;
    nx=length(x);

    % define BaS and its initial and final states
    w=0;
    for ii=1:length(h)
        w=w+1/h{ii}(x);
    end
    w = matlabFunction(w,'Vars',{x});
    w0=w(x0);
    wf=w(xf);
    
    % define safety equations at the next state h(xk+1)=h(f(xk,uk)):
    hplus=cell(length(h));
    for ii=1:length(hplus)
        hplus{ii}= @(x,u) h{ii}(xplus(x,u));
    end
    % find gradient of h
        hx=jacobian(h{1}(x),x);

%     hx= matlabFunction(h,'Vars',{x});
    beta_prime = -1/(h{1}(x))^2;
    % DBaS state equation
    x = sym('x',[nx+n_bas 1]); %n_bas is number of barrier states;
    w_d=0;
    for ii=1:length(hplus)
        w_plus= x(nx+1) + dt*(beta_prime*hx*xdot);
    end
    
    f_w= w_plus;
        
    fx_w=jacobian(f_w,x);
    fu_w=jacobian(f_w,u); 

    % symbolic to function handle
    f_w= matlabFunction(f_w,'Vars',{x,u});
    fx_w= matlabFunction(fx_w,'Vars',{x,u});
    fu_w= matlabFunction(fu_w,'Vars',{x,u});
end

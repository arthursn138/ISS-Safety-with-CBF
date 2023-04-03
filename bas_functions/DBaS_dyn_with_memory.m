%% Dynamics of Discrete Barrier States
function [f_w,fx_w,fu_w,w0,wf]=DBaS_dyn_with_memory(x,u,x0,xf,f_dyn,h,alpha,n_bas)
    xplus=@(x,u) f_dyn(x,u);
    nx=length(x);

    % define DBaS and its initial and final states
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
    % DBaS state equation
    w_plus=0;
    for ii=1:length(hplus)
        w_plus=w_plus+1/hplus{ii}(x,u);
    end
    
    x = sym('x',[nx+n_bas 1]); %n_bas is number of barrier states;
    % w = x(nx+1), the new state
    
    f_w= alpha*x(nx+1) + (1-alpha)*w_plus;

    fx_w=jacobian(f_w,[x]);
    fu_w=jacobian(f_w,u); 

    % symbolic to function handle
    f_w= matlabFunction(f_w,'Vars',{x,u});
    fx_w= matlabFunction(fx_w,'Vars',{x,u});
    fu_w= matlabFunction(fu_w,'Vars',{x,u});
end

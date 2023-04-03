%% Dynamics of Discrete Barrier States
function [f_w,fx_w,fu_w,w0,wf,fxx_w,fxu_w,fuu_w]=DBaS_dyn(x,u,x0,xf,f_dyn,h,ddp_2nd_order)
    xplus=@(x,u) f_dyn.f(x,u);
    % define DBaS and its initial and final states
    w=0;
    for ii=1:length(h)
        w=w+1/h{ii}(x);
%         w=w+1/h{ii}(x) - 1/h{ii}(xf);
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
%         w_plus=w_plus+1/(1/w(x) + hplus{ii}(x,u) - h{ii}(x));
    end
    
    f_w= w_plus;
    syms xplus % w = xplus, the new state
    xbar = [x;xplus];
    fx_w=jacobian(f_w,xbar);
    fu_w=jacobian(f_w,u); 
%     for ii=1:length(xbar)
%     fxx_w=[]; fxu_w=[];
    fxx_w=jacobian(fx_w, xbar);
    fxu_w=jacobian(fu_w,xbar);
    fuu_w=jacobian(fu_w,u);
%     end
    fxu_w=fxu_w.';
    % symbolic to function handle
    f_w= matlabFunction(f_w,'Vars',{x,u});
    fx_w= matlabFunction(fx_w,'Vars',{x,u});
    fu_w= matlabFunction(fu_w,'Vars',{x,u});
    if ddp_2nd_order
        fxx_w= matlabFunction(fxx_w,'Vars',{x,u});
        fxu_w= matlabFunction(fxu_w,'Vars',{x,u});
        fuu_w= matlabFunction(fuu_w,'Vars',{x,u});
    end
end

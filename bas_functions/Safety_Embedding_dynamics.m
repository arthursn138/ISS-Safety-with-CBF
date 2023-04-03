function [fbar,fbarx,fbaru,fbarxx,fbarxu,fbaruu,xbar0,xbarf,nbar]=Safety_Embedding_dynamics(x,u,x0,xf,f_dyn,w0,wf,f_w,fx_w,fu_w,fxx_w,fxu_w,fuu_w,deriv_bool,ddp_2nd_order)
    n=length(x); m=length(u);
    % augment bas dynamics
    fbar=@(x,u)[f_dyn.f(x,u);f_w(x,u)];
    % augment gradients
    [fbarx, fbaru,fbarxx,fbarxu,fbaruu]=deal([]);
    if deriv_bool
        % dynamics gradients
        fbarx=@(x,u)[f_dyn.fx(x,u) zeros(n,1);fx_w(x,u)];
        fbaru=@(x,u) [f_dyn.fu(x,u);fu_w(x,u)];
        if ddp_2nd_order
%             fbarxx=f_dyn.fxx;
%             fbarxx{n+1}=fxx_w;
            syms xplus; x=[x;xplus]; fx_temp=fbarx(x,u);
            for ii=1:n+1
                fxx_temp = jacobian(fx_temp(ii,:), x);
                fbarxx{ii}=matlabFunction(fxx_temp,'Vars',{x,u});
                fxu_temp = jacobian(fx_temp(ii,:), u).';
                fbarxu{ii}=matlabFunction(fxu_temp,'Vars',{x,u});
            end
        fu_temp = fbaru(x,u);
        for ii=1:m
                fuu_temp = jacobian(fu_temp(ii,:), u);
                fbaruu{ii}=matlabFunction(fuu_temp,'Vars',{x,u});
        end
        end
    end
    
    xbar0=[x0;w0]; xbarf=[xf;wf];
    nbar=length(xbar0);
end
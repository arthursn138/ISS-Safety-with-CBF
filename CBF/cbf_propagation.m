% CBF propagation
function [dx u]=cbf_propagation(t,x,f,g,P,alpha_clf,cbf_x,alpha_cbf,p,u0,umin,umax,options)
    %Dynamics
    f_clf=f(x);
    g_clf=g(x);
    cbf_x_val=cbf_x{1}(x);
    alpha_cbf=alpha_cbf(x);
    
    [u_clf]=fmincon(@(u_clf) cbf_clf_objective(u_clf,p),u0,[],[],[],[],umin,umax,@(u_clf) cbf_clf_constraint(u_clf,f_clf,g_clf,x,P,alpha_clf,cbf_x_val,alpha_cbf),options);
    
    dx=f_clf+g_clf*u_clf(1);
end

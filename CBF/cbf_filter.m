function [u]=cbf_filter(u,f,g,hx,alpha)
% [T_zcbf_clf,X_zcbf_clf] = ode45(@(t,x) cbf_clf_ode(t,x,f,g,P,alpha_clf,hx,alpha_cbf,relaxation_penalty,u0_clf_cbf,clf_cbf_min,clf_cbf_max,options), tspan, x0);

% composite barrier:
options = optimoptions('fmincon','Display','off');

[u]=fmincon(@(u) cbf_objective(u),u0,[],[],[],[],[],[],@(u) cbf_constraint(u,f,g,hx,alpha),options);


end

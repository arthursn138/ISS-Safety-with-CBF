function [l,l_x, l_u, l_xu, l_xx, l_uu] = run_quad_cost_penalty(x,u,Q,R,xf,beta,beta_x,beta_xx,Q_bf,deriv_bool)
% running quadratic (in x and u) cost    
    e= (x-xf);

    l = 0.5*( e'*Q*e + u.'*R*u + Q_bf*(beta(e))^2);
    if deriv_bool
        % cost derivatives (grad and Hess)
        l_x=Q*e + Q_bf*beta_x(e).';
        l_u=R*u;
        l_xu=zeros(length(x),length(u));
        l_xx=Q + Q_bf*beta_xx(e);
        l_uu=R;
    else
        [l_x, l_u,l_xu, l_xx, l_uu]=deal([]);
    end
end
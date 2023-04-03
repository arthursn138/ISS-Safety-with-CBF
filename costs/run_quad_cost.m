function [l,l_x, l_u, l_xu, l_xx, l_uu] = run_quad_cost(x,u,Q,R,xf,deriv_bool)
% running quadratic (in x and u) cost    
    l = 0.5*( (x-xf)'*Q*(x-xf) + u.'*R*u );
    if deriv_bool
        % cost derivatives (grad and Hess)
        l_x=Q*(x-xf);
        l_u=R*u;
        l_xu=zeros(length(x),length(u));
        l_xx=Q;
        l_uu=R;
    else
        [l_x, l_u,l_xu, l_xx, l_uu]=deal([]);
    end
end
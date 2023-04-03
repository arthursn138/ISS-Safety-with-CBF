function [l_terminal,l_term_x,l_term_xx] = terminal_quad_cost(x,xf,S,deriv_bool)
% running quadratic (in x and u) cost    
    l_terminal = 0.5*(x-xf)'*S*(x-xf);
    if deriv_bool  % no gradients
        l_term_x = S*(x - xf);
        l_term_xx = S;
    else
        [l_term_x, l_term_xx] = deal([]);
    end
end
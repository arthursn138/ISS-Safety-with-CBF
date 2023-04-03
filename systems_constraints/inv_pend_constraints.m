function [h, vel_limit, n_bas]=inv_pend_constraints()
vel_limit=5;
h{1} = @(x) vel_limit^2 - x(2)^2;
n_bas=1;
end
    

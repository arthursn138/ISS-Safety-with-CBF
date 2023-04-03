function [h, control_limit]=control_constraints()
control_limit=15;
h{1} = @(x,u) control_limit^2 - u^2;
end
    

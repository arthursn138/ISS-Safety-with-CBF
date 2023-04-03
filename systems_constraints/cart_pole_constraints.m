function [h, cart_limit]=cart_pole_constraints()
cart_limit=1.5;
h{1} = @(x) cart_limit^2 - x(1)^2;
end
    

 % CBF constraints
function [c,ceq] = cbf_constraint(u,f,g,hx,alpha)
c = [-(hx*(f+g*u)+alpha)];
ceq = [];
end
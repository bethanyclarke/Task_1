function [A] = rcross(x)
% Given a column vector x, produces the matrix A such that A*v = cross(v,x);

A = [0,x(3),-x(2);-x(3),0,x(1);x(2),-x(1),0];


end


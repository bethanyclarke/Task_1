function [out] = dexpinv(v,u)
% dexpinv equation for (R^3, \times).

theta = norm(v);

if theta < 10^-6
    fac = -1/12;
else
    fac = (0.5*theta*cot(0.5*theta) - 1)/(theta^2);
end

out = u - 0.5*cross(v,u) - fac*cross(v,cross(v,u));

end


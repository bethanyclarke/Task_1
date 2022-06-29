function [p] = unit_q_sqrt(q)

if q(1)<-0.999999 % Avoid the p = -1 case.
    p = [1,0,0,0];
else
    a = sqrt(0.5*(q(1)+1));
    p = [a,[q(2),q(3),q(4)]/(2*a)];
end

end


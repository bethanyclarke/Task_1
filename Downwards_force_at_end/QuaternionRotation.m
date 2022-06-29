function [out] = QuaternionRotation(q,v)

out = v(1)*[1 - 2*(q(3)^2 + q(4)^2); 2*(q(2)*q(3) + q(1)*q(4)); 2*(q(2)*q(4) - q(1)*q(3))] ...
                + v(2)*[2*(q(2)*q(3) - q(1)*q(4));1 - 2*(q(2)^2 + q(4)^2);2*(q(3)*q(4) + q(1)*q(2))] ...
                + v(3)*[2*(q(2)*q(4) + q(1)*q(3));2*(q(4)*q(3) - q(1)*q(2));1 - 2*(q(3)^2 + q(2)^2)];
end

function [C] = J_block_C(Fibre,fac)

Np = Fibre.Np;

Q = Fibre.Q;

C = zeros(3*(Np-1));

C(1:3,1:3) = fac*rcross(QuaternionRotation(Q(2,1:4),[1;0;0]));

for k=3:Np
    
    i = 3*(k-2);
    
    C(i+1:i+3,1:i) = C(i-2:i,1:i);
    
    C(i+1:i+3,i-2:i+3) = C(i+1:i+3,i-2:i+3) + fac*[rcross(QuaternionRotation(Q(k-1,1:4),[1;0;0])),rcross(QuaternionRotation(Q(k,1:4),[1;0;0]))];
    
end

end


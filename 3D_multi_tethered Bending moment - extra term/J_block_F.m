function [F] = J_block_F(Fibre,dt,mu)

Np = Fibre.Np;

dL = Fibre.DeltaL;

R = Fibre.R;

U = Fibre.U;

Q = Fibre.Q;

F = zeros(3*Np,3*(Np-1));

for i=1:Np-1
    
    pos = 3*(i-1);
    
    fac = dt*dL/(24*pi*mu*R^3);
    
    tmat = rcross(QuaternionRotation(Q(i,1:4),[1;0;0]))';
    umat = rcross(U(1:3,i))';
    
    F(pos+1:pos+3,pos+1:pos+3) = fac*(tmat - 0.5*umat*tmat + umat^2*tmat/12);
    
    fac = dt*dL/(24*pi*mu*R^3);
    
    tmat = rcross(QuaternionRotation(Q(i+1,1:4),[1;0;0]))';
    umat = rcross(U(1:3,i+1))';
    
    F(pos+4:pos+6,pos+1:pos+3) = fac*(tmat - 0.5*umat*tmat + umat^2*tmat/12);
    
end

F = F(4:end,:);

end


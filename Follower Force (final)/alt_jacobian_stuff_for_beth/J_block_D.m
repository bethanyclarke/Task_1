function [D] = J_block_D(Fibre,mu,dt)

Np = Fibre.Np;

fac = zeros(1,3*(Np-1));

heights = Fibre.X(3,:)/Fibre.R;

for n=2:Np
    
    id = 3*(n-2);
    
    h = heights(n);
    
    fac(id+1) = 1 - (9/h -2*h^(-3) + h^(-5))/16;
    fac(id+2) = fac(id+1);
    fac(id+3) = 1 - (9/h -4*h^(-3) + h^(-5))/8;
    
end

fac = -2*dt*fac/(18*pi*mu*Fibre.R);

D = diag(fac);

for n=2:Np-1
    
    id = 3*(n-2);
    
    D(id+1,id+4) = -D(id+1,id+1);
    D(id+2,id+5) = -D(id+2,id+2);
    D(id+3,id+6) = -D(id+3,id+3);
    
end

end


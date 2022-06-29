function [check,errorvec] = ConstraintCheck(Filaments,dt,omega)

% This code checks whether or not the constraints we impose on the system
% are satisfied.

check = 0;

TOL = 10^(-4);

Nf = length(Filaments);

c1 = 4/3;
c2 = 1/3;
c3 = 2/3;

for i=Nf:-1:1
    
    N(1,i) = 6*Filaments(i).Np;
    
end

N = [0,cumsum(N)];

errorvec = zeros(N(end),1);

for n=1:Nf
    
    dL = Filaments(n).DeltaL;
    
    V = Filaments(n).V;
    
    Q = Filaments(n).Q;
    
    U = Filaments(n).U;
    
    Vtarg = V(1:3,1);
    
    for i=1:Filaments(n).Np
        
        if i==1
            
            pos_error = Vtarg;
            U_error = V(4:6,1) - [0;0;omega];
            
        else
            
            tangent_bdf = zeros(3,1);
            
            for pid=[i-1,i]
                
                oldQ = QuaternionProduct(qexp(U(4:6,pid)),Q(pid,5:8)); % Undo the pervious orientation update.
                
                tangent_bdf = tangent_bdf + QuaternionRotation(Q(pid,1:4),[1;0;0]) ...
                    - c1*QuaternionRotation(Q(pid,5:8),[1;0;0]) ...
                    + c2*QuaternionRotation(oldQ,[1;0;0]);
                
            end
            
            Vtarg = Vtarg + ((3*dL)/(4*dt))*tangent_bdf;
            
            pos_error = V(1:3,i) - Vtarg;
            
            U_error = U(1:3,i) - (1/3)*U(4:6,i) - ...
            c3*dt*dexpinv(U(1:3,i),V(4:6,i));
            
        end
        
        loc = N(n) + 3*(i-1);
        errorvec(loc+1) = pos_error(1);
        errorvec(loc+2) = pos_error(2);
        errorvec(loc+3) = pos_error(3);
        errorvec(loc+1+3*Filaments(n).Np) = U_error(1);
        errorvec(loc+2+3*Filaments(n).Np) = U_error(2);
        errorvec(loc+3+3*Filaments(n).Np) = U_error(3);
        
        e = max(abs([pos_error;U_error]));
        
        if e > TOL
            
            check = 1;
            
        end
        
    end
    
end

end


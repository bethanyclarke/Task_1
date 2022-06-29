function [check,errorvec] = ConstraintCheck(Filaments,dt,nt)

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

if nt==1 % Backwards Euler
    
    for n=1:Nf
        
        dL = Filaments(n).DeltaL;
        
        V = Filaments(n).V;
        
        Q = Filaments(n).Q;
        
        U = Filaments(n).U;
        
        Xcurr = Filaments(n).X(:,1);
        
        Xm1 = Filaments(n).Xm1;
        
        for i=1:Filaments(n).Np
            
            if i==1
                
                pos_error = V(1:3,i);
                U_error = U(1:3,i) - dt*dexpinv(U(1:3,i),V(4:6,i));
                
            else
                
                for pid=[i-1,i]

                    Xcurr = Xcurr + 0.5*dL*QuaternionRotation(Q(pid,1:4),[1;0;0]);
                    
                end
                
                pos_error = Xcurr - Xm1(:,i) - dt*V(1:3,i);
                
                U_error = U(1:3,i) - dt*dexpinv(U(1:3,i),V(4:6,i));
                
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
    
else % BDF2
    
    for n=1:Nf
        
        dL = Filaments(n).DeltaL;
        
        V = Filaments(n).V;
        
        Q = Filaments(n).Q;
        
        U = Filaments(n).U;
        
        Xcurr = Filaments(n).X(:,1);
        
        Xm1 = Filaments(n).Xm1;
        
        Xm2 = Filaments(n).Xm2;
        
        for i=1:Filaments(n).Np
            
            if i==1
                
                pos_error = V(1:3,i);
                U_error = U(1:3,i) - (1/3)*U(4:6,i) - ...
                    c3*dt*dexpinv(U(1:3,i),V(4:6,i));
                
            else
                
                for pid=[i-1,i]
                    
                    Xcurr = Xcurr + 0.5*dL*QuaternionRotation(Q(pid,1:4),[1;0;0]);
                    
                end

                pos_error = Xcurr - c1*Xm1(:,i) + c2*Xm2(:,i) - c3*dt*V(1:3,i);
                
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

end


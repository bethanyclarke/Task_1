function [] = StericForces(Filaments)

Nf = length(Filaments);

for i=1:Nf
    
    Xi = Filaments(i).X;
    
    Ri = Filaments(i).R;
    
    Fi = Filaments(i).F;
    
    for j=i:Nf
        
        Xj = Filaments(j).X;
        
        Rj = Filaments(j).R;
        
        Fj = Filaments(j).F;
        
        x = norm(Filaments(i).CentreOfMass - Filaments(j).CentreOfMass);
        d = 0.5*(Filaments(i).Np * Filaments(i).DeltaL + Filaments(j).Np * Filaments(j).DeltaL);
        
        if x <= d
            
            for m=1:Filaments(i).Np
                
                Xm = Xi(:,m);
                
                r = Ri;
                
                for n=1:Filaments(j).Np
                    
                    diff = Xm - Xj(:,n);
                    
                    dist2 = diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3);
                    
                    r2 = (r + Rj)^2;
                    
                    chi = 1.21 * r2;
                    
                    if dist2 < chi
                        
                        fac = 15*pi*((chi - dist2)/(chi - r2))^4;
                        
                        Force = fac*diff;
                        
                        Fi(1:3,m) = Fi(1:3,m) + Force;
                        
                        if i~=j % So that we don't cancel out barrier forces in the same filament.
                            
                            Fj(1:3,n) = Fj(1:3,n) - Force;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        Filaments(j).F = Fj;
        
    end
    
    Filaments(i).F = Fi;
    
end

end


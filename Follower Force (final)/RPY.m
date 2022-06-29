function [] = RPY(Filaments,mu)

N = length(Filaments);

Nparticles = 0;

for i=1:N
    
    Nparticles = Nparticles + Filaments(i).Np;
    
end

X = zeros(1,3*Nparticles);
Forces = zeros(1,3*Nparticles);
Torques = zeros(1,3*Nparticles);
R = zeros(1,Nparticles);

myN = 0;

for i=1:N
    
    X(3*myN+1:3*myN+3*Filaments(i).Np) = Filaments(i).X(1:end);
    
    temp = Filaments(i).F(1:3,:);
    Forces(3*myN+1:3*myN+3*Filaments(i).Np) = temp(1:end);
    
    temp = Filaments(i).F(4:6,:);
    Torques(3*myN+1:3*myN+3*Filaments(i).Np) = temp(1:end);
    
    R(myN+1:myN+Filaments(i).Np) = Filaments(i).R*ones(1,Filaments(i).Np);
    
    myN = myN + Filaments(i).Np;
    
end

% Pass to the mex function.

% [V,Omega] = mexRPY(X,Forces,Torques,R,mu);
[V,Omega] = mexRPYwall(X,Forces,Torques,R(1),mu); % All radii must be the same in this case

% Place the velocities into the filament objects.

myN = 0;

for i=1:N
    
    Filaments(i).V = [reshape(V(3*myN+1:3*myN+3*Filaments(i).Np),3,Filaments(i).Np);reshape(Omega(3*myN+1:3*myN+3*Filaments(i).Np),3,Filaments(i).Np)];
    
    myN = myN + Filaments(i).Np;
    
end

end % End function.


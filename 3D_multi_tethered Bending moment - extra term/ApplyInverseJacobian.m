function [out] = ApplyInverseJacobian(v,Filaments,Cmat,Dmat,iter)

for i=length(Filaments):-1:1
    
    N(i) = 6*Filaments(i).Np;
    
end

N = [0,cumsum(N)];

out = zeros(N(end),1);

for i=1:length(Filaments)
    
    out(N(i)+1:N(i+1)) = Filaments(i).InvertLocalBlock(v(N(i)+1:N(i+1)));
    
end

for i=1:iter
    
    out = out + (Dmat(:,i)' * v)*Cmat(:,i);
    
end

end


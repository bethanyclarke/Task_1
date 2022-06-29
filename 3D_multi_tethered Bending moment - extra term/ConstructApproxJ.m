function [J] = ConstructApproxJ(Fibre,dt,mu)

Np = Fibre.Np; % The number of particles in the filament.

J = zeros(6*Np);

a = Fibre.R;

dL = Fibre.DeltaL;

% We start with the entries of the Jacobian corresponding to the update
% equation for the position of the first particle.

j = 3*Np + 3;
fac = 1/(6*pi*mu*a);

J(1,1) = fac; J(2,2) = fac; J(3,3) = fac;
J(1,j+1) = -fac; J(2,j+2) = -fac; J(3,j+3) = -fac;

% Next, we provide the terms related to the derivative of the velocity
% constraints with respect to the Lie algebra elements.

J(4:3*Np,7:3*Np+3) = J_block_C(Fibre,-0.75*Fibre.DeltaL/dt);

% Now, we produce the block relating the velocity constraints to the
% Lagrange multipliers.

J(4:3*Np,3*Np+4:end) = J_block_D(Fibre,mu);

for m=1:Np-1
    J(3*m+1:3*m+3,1:3) = -eye(3)/(6*pi*a*mu);
end

% Next, the block encoding the dependence of the Lie algebra update
% equations on the Lie algebra elements.

J(3*Np+4:end,7:3*Np+3) = J_block_E(Fibre,dt,mu);

% Finally, we produce the derivatives of the Lie algebra update equations
% with repsect to the Lagrange multipliers.

J(3*Np+4:end,3*Np+4:end) = J_block_F(Fibre,dt,mu);

J(3*Np+1:3*Np+3,3*Np+4:3*Np+6) = -dL*rcross(QuaternionRotation(Fibre.Q(1,1:4),[1;0;0]))'/(16*pi*mu*a^3);

J(3*Np+1:3*Np+3,4:6) = eye(3)/(8*pi*mu*a^3);

end


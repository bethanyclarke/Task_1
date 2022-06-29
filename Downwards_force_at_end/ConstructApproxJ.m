function [J] = ConstructApproxJ(Fibre,dt,mu)

Np = Fibre.Np; % The number of particles in the filament.

J = zeros(6*Np);

a = Fibre.R;

dL = Fibre.DeltaL;

% We start with the entries of the Jacobian corresponding to the update
% equation for the position of the first particle.

j = 3*Np + 3;
h = Fibre.X(3,1)/a; % normalised height of the first particle above the wall.
fac1 = (1 - (9/h -2*h^(-3) + h^(-5))/16)/(6*pi*mu*a);
fac2 = (1 - (9/h -4*h^(-3) + h^(-5))/8)/(6*pi*mu*a);

J(1,1) = fac1; J(2,2) = fac1; J(3,3) = fac2;
J(1,j+1) = -fac1; J(2,j+2) = -fac1; J(3,j+3) = -fac2;

% Next, we provide the terms related to the derivative of the position
% constraints with respect to the Lie algebra elements.

J(4:3*Np,4:3*Np+3) = J_block_C(Fibre,0.5*dL);

% Now, we produce the block relating the velocity constraints to the
% Lagrange multipliers.

J(4:3*Np,3*Np+4:end) = J_block_D(Fibre,mu,dt);

% Next, the block encoding the dependence of the Lie algebra update
% equations on the Lie algebra elements.

J(3*Np+1:end,4:3*Np+3) = J_block_E(Fibre,dt,mu);

% Finally, we produce the derivatives of the Lie algebra update equations
% with repsect to the Lagrange multipliers.

J(3*Np+1:end,3*Np+4:end) = J_block_F(Fibre,dt,mu);

end


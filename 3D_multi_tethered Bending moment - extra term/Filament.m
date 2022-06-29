classdef Filament < handle
    
    % N.B. This is a handle class and does NOT have a copy constructor.
    % Assignment to a Filament type will cause both variables to be a
    % reference to the same memory location.
    
    properties
        
        Np; % Number of beads comprising the filament.
        
        Kb; % The elastic bending modulus of the filament.
        
        Kt; % The elastic twisting modulus of the filament.
        
        DeltaL; % The centreline distance between adjacent beads.
        
        StrainTwist; % The strain-twist vector encoding the intrinsic curvature and twist of the filament.
        
        Bu; % The Buckling number.
        
        omega; % Rotation rate.
        
        R; % The common radius of the constituent particles.
        
        X; % Positions of the particles comprising the filament.
        
        Q; % Orientation quaternions for the particles comprising the filament.
        
        U; % Lie algebra elements associated with the orientation of the constituent particles.
        
        V; % Velocities (both angular and translational) of the particles.
        
        F; % Forces and torques on the particles.
        
        Lambda; % The collection of Lagrange multipliers associated with the inextensibility of the filament.
        
        TetherLam; % The vector Lagrange multiplier associated with the tethering constraint.
        
        RotLam; % The vector Lagrange multiplier associated with the imposed rotation.
        
        Lmat; % The lower-triangular part of the Jacobian L-U decomposition.
        
        Umat; % The upper-triangular part of the Jacobian L-U decomposition.
        
    end
    
    methods
        
        function obj = Filament(varargin)
            
            if nargin==0
                
                % Empty constructor doesn't actually need to do anything,
                % it's only so that we can pre-allocate arrays to contain
                % Filaments.
                
            else
                
                N = varargin{1};
                
                obj.Np = N;
                
                obj.R = 1/(2.2*N);
                
                obj.X = zeros(3,N);
                
                obj.Q = zeros(N,8);
                
                obj.U = zeros(9,N);
                
                obj.V = zeros(6,N);
                
                obj.F = zeros(6,N);
                
                obj.DeltaL = 2.2*obj.R;
                
                obj.Lambda = zeros(3,N-1);
                
                obj.TetherLam = zeros(3,1);
                
                obj.RotLam = zeros(3,1);
                
            end
            
        end
        
        function RobotArm(obj)
            
            Xtemp = obj.X;
            
            Qtemp = obj.Q;
            
            dL = obj.DeltaL;
            
            for i=2:obj.Np
                
                Xtemp(:,i) = Xtemp(:,i-1) + ...
                    0.5*dL*(QuaternionRotation(Qtemp(i-1,1:4),[1;0;0]) + ...
                    QuaternionRotation(Qtemp(i,1:4),[1;0;0]));
                
            end
            
            obj.X = Xtemp;
            
        end
        
        function InitialSetup(obj,FirstBeadPos,StrainTwist,Bu,omega,mu,Nf,h,dt)
            
            obj.X(:,1) = FirstBeadPos; % Because the filaments are tethered,
            % the position of the first bead is set here and never changed again.
            
            L = obj.Np * obj.DeltaL;
            
            a = obj.R;
            
            fac = 2^-0.5;
            
            obj.Q = [fac*ones(obj.Np,1),zeros(obj.Np,1),-fac*ones(obj.Np,1),zeros(obj.Np,1),...
                fac*ones(obj.Np,1),zeros(obj.Np,1),-fac*ones(obj.Np,1),zeros(obj.Np,1)];
            % The orientation quaternions are initialised such that the filaments point up the z-axis by default.
            
            obj.RobotArm;
            
            obj.StrainTwist = StrainTwist;
            
            obj.Bu = Bu;
            
            obj.omega = omega;
            
            perp = 4*pi*mu/log(L/a);
            
            obj.Kb = perp*omega*L^4*(a/h)^2*(Nf-1)/obj.Bu;
            
            obj.Kt = obj.Kb/10;
            
            obj.U(1:3,1) = [0;0;0];
            
            obj.U(4:6,1) = [0;0;0];
            
            obj.U(7:9,1) = [0;0;0];
            
        end
        
        function InitialGuess(obj)
            
            Utemp = obj.U;
            
            Qtemp = obj.Q;
            
            for i=1:obj.Np
                
                % Guess Lie algebra elements, and hence quaternions, of all
                % particles.
                
                if i>1

                    Utemp(1:3,i) = 2*Utemp(4:6,i) - Utemp(7:9,i);
                    
                end
                
                Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),Qtemp(i,5:8));
                
            end
            
            obj.U = Utemp;
            
            obj.Q = Qtemp;
            
        end
        
        function InternalForcesAndTorques(obj)
            
            N = obj.Np;
            
            ForcesAndTorques = zeros(6,N);
            
            Qtemp = obj.Q;
            
            Lam = obj.Lambda;
            
            dL = obj.DeltaL;
            
            ST = obj.StrainTwist;
            
            Bend = obj.Kb;
            
            Twist = obj.Kt;
            
            Tether = obj.TetherLam;
            
            Rot = obj.RotLam;
            
            for i=1:obj.Np-1
                
                % Constraint forces and torques
                
                ForcesAndTorques(1:3,i) = ForcesAndTorques(1:3,i) - Lam(:,i);
                
                ForcesAndTorques(1:3,i+1) = ForcesAndTorques(1:3,i+1) + Lam(:,i);
                
                t = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);
                
                ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) - ...
                    0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) - t(1)*Lam(3,i);t(1)*Lam(2,i) - t(2)*Lam(1,i)];
                
                t = QuaternionRotation(Qtemp(i+1,1:4),[1;0;0]);
                
                ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - ...
                    0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) - t(1)*Lam(3,i);t(1)*Lam(2,i) - t(2)*Lam(1,i)];
                
                % Elastic torques
                
                q = MidpointQ(Qtemp(i,1:4),Qtemp(i+1,1:4));
                
                dqds = (Qtemp(i+1,1:4) - Qtemp(i,1:4))/dL;
                
                b = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
                b = 2*b(2:4);
                
                M = QuaternionRotation(q,[Twist * (b(1) - ST(3)); ...
                    Bend * (b(2) - ST(1)); Bend * (b(3) - ST(2))]);
                
                ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) + M;
                
                ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - M;
                
            end
            
            % Tethering force and torque
            
            u = obj.U(1:3,1);
            
            ForcesAndTorques(1:3,1) = ForcesAndTorques(1:3,1) + Tether;
            
            theta = (u(1)*u(1) + u(2)*u(2) + u(3)*u(3))^0.5;
            
            if theta<10^-6
                fac = -1/12;
            else
                fac = (0.5*theta*cot(0.5*theta) - 1)/(theta^2);
            end
            
            TetherTorque = Rot - 0.5*cross(Rot,u) - fac*cross(cross(Rot,u),u);
            
            ForcesAndTorques(4:6,1) = ForcesAndTorques(4:6,1) + TetherTorque;
            
            obj.F = ForcesAndTorques;
            
        end
        
        function com = CentreOfMass(obj)
            
            com = mean(obj.X,2);
            
        end
        
        function DecomposeJacobian(obj,dt,mu)
            
            [obj.Lmat,obj.Umat] = lu(ConstructApproxJ(obj,dt,mu));
            
        end
        
        function out = InvertLocalBlock(obj,v)
            
            out = obj.Lmat\v;
            
            out = obj.Umat\out;
            
        end
        
        function ApplyUpdate(obj,u)
            
            N = obj.Np;
            
            obj.TetherLam = obj.TetherLam + u(1:3);
            
            obj.RotLam = obj.RotLam + u(4:6);
            
            Utemp = obj.U;
            
            Qtemp = obj.Q;
            
            Lam = obj.Lambda;
            
            for i=2:N
                
                Utemp(1:3,i) = Utemp(1:3,i) + u(3*i+1:3*(i+1));
                
                Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),Qtemp(i,5:8));
                
            end
            
            for i=1:N-1
                
                Lam(:,i) = Lam(:,i) + u(3*(N+i)+1:3*(N+i+1));
                
            end
            
            obj.U = Utemp;
            
            obj.Q = Qtemp;
            
            obj.Lambda = Lam;
            
        end
        
        function EndOfStepUpdate(obj)
            
            Utemp = obj.U;
            
            Qtemp = obj.Q;
            
            Utemp(7:9,:) = Utemp(4:6,:);
            
            Utemp(4:6,:) = Utemp(1:3,:);
            
            for i=1:obj.Np
                
                Qtemp(i,5:8) = QuaternionProduct(qexp(Utemp(1:3,i)),Qtemp(i,5:8));
                
            end
            
            obj.U = Utemp;
            
            obj.Q = Qtemp;
            
        end
        
    end
    
end


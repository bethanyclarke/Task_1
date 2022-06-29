function [] = main(Bu)

Nf = 2; % Number of filaments.
Np = 50; % Number of particles comprising each filament.

mu = 1; % Fluid viscosity.

omega = 0.1; % Rotation rate of the filaments.

h = 0.1; % Diameter of the circle along the edge of which the filaments are tethered.
% We should have 1/(2.2*Np) << h << 1.

theta = 2*pi/Nf;

StepsPerPeriod = 50;
TotalSteps = 400*StepsPerPeriod;
dt = 2*pi/(StepsPerPeriod*omega);

VideoName = sprintf('Bu_equals_%i_video.avi',Bu);
video = VideoWriter(VideoName);
open(video);

for i=Nf:-1:1
    
    Filaments(i) = Filament(Np);
    Filaments(i).InitialSetup([0.5*h*cos((i-1)*theta);0.5*h*sin((i-1)*theta);0],[0,0,0],Bu,omega,mu,Nf,h,dt);
    
    N(i) = 6*Filaments(i).Np;
    
end

Tb = 4*pi*mu/(log(2.2*Np) * Filaments(1).Kb);

N = [0,cumsum(N)];

Nbroy = N(end);

for Steps=1:TotalSteps
    
    Cmat = zeros(Nbroy,20); Dmat = zeros(Nbroy,20);
    
    for i=1:Nf
        
        Filaments(i).InitialGuess;
        Filaments(i).RobotArm;
        Filaments(i).InternalForcesAndTorques;
        
    end
    
    StericForces(Filaments);
    
    RPY(Filaments,mu);
    
    [Check,Error] = ConstraintCheck(Filaments,dt,omega);

    for i=1:Nf
        
        Filaments(i).DecomposeJacobian(dt,mu);
        
    end
    
    BroydenIter = 0;
    
    while Check==1
        
        Update = ApplyInverseJacobian(-Error,Filaments,Cmat,Dmat,BroydenIter);

        for i=1:Nf
            
            Filaments(i).ApplyUpdate(Update(N(i)+1:N(i+1)));
            Filaments(i).RobotArm;
            Filaments(i).InternalForcesAndTorques;
            
        end
        
        StericForces(Filaments);
        
        RPY(Filaments,mu);
        
        [Check,NewError] = ConstraintCheck(Filaments,dt,omega);
        
        Svec = ApplyInverseJacobian(-NewError,Filaments,Cmat,Dmat,BroydenIter);
        
        Diff = NewError - Error;
        
        BroydenIter = BroydenIter + 1;
        
        fac = norm(Diff);
        
        Cmat(:,BroydenIter) = Svec/fac;
        
        Dmat(:,BroydenIter) = Diff/fac;
        
        Error = NewError;

    end
    
    for i=1:Nf
        
        Filaments(i).EndOfStepUpdate;
        
    end
    
    fprintf('Step %i/%i required %i Broyden iterations.\n',Steps,TotalSteps,BroydenIter);
    
    if mod(Steps,20)==1
        
        % Plotting
        
        [s1,s2,s3] = sphere;
        
        for i=1:Nf
            
            MyColour = ((i-1)/Nf)*[1,1,1];
            
            a = Filaments(i).R;
            
            x = Filaments(i).X;
            
            for j=1:Filaments(i).Np
                
                surf(x(1,j)+a*s1,x(2,j)+a*s2,x(3,j)+a*s3,'FaceColor',MyColour);
                hold on;
                
            end
 
        end
       
        axis equal;
        view(0,0);
        
        PlotName = sprintf('Bu = %i, t/Tb = %.5g',Bu,Steps*dt/Tb);
        title(PlotName);
        
        CurrFrame = getframe(gcf);
        writeVideo(video,CurrFrame);
        
        pause(0.1);
        hold off;
        
    end
    
end

close(video);

end


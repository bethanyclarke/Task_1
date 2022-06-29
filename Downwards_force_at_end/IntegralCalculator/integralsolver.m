% Solving integral, using fzero (less good then Tims way!). Ended up using untitled.m instead
N = 20; % Number of particles
dL = 2.2; % Distance between segements
L = N * dL - 2; % Length of filament
Bend = 1800; % Bending modulus
F = 1.93*Bend/L^2; % Horizontal force



% D = zeros(1, 2);
% 
% L = 44;
% f_fac = 1/sqrt(500); % 0.5*KB/f = 0.5*KB/(0.001*KB) = 500.
% 
% integral = @(theta0, theta) abs(real((2/sqrt(cos(theta0)-1))*(ellipticF(pi/4, -2/(cos(theta0)-1)) - ellipticF(0.5*theta, -2/(cos(theta0)-1)))));
% 
% D(:,1) = 1./NSEG';
% 
% figure;
% hold on;
% 
% for n=1:length(NSEG)
%        
%     if exist(theta_file_name, 'file')
%                 
%     else
%         
%         fprintf('No theta data file was found for a %i-segment filament. Calculating analytical values...\n', NSEG(n));
%         
%         theta0 = 8.054788116455077e-01; % Previously calculated solution to integral(theta0, theta0) = L/sqrt(0.5*KB/f) = 44/sqrt(500).
%         
%         theta(NSEG(n)) = theta0; % theta(k) = \theta(s = (k-1)*L/(NSEG(n)-1)).
%         theta(1) = pi/2;
%         
%         for k=2:NSEG(n)-1
%             
%             upper_bound = theta(k-1);
%             lower_bound = theta0;
%             
%             theta(k) = 0.5*(upper_bound + lower_bound);
%             
%             s = (k-1)*L/(NSEG(n)-1);
%             
%             error = s*f_fac - integral(theta0, theta(k));
%             
%             while abs(error) > 1e-9
%                 
%                 if error > 0
%                     
%                     upper_bound = theta(k);
%                     
%                 else
%                     
%                     lower_bound = theta(k);
%                     
%                 end
%                 
%                 theta(k) = 0.5*(upper_bound + lower_bound);
%                 
%                 error = s*f_fac - integral(theta0, theta(k));
%                 
%             end
%             
%             fprintf('Calculated theta(%i) = %g\n', k, theta(k));
%             
%         end
%         
%         dlmwrite(theta_file_name, theta, 'precision', '%.12e');
%         
%     end
%     
%     filD = load(sprintf('%isegs.dat', NSEG(n)));
%     
%     plot((0 : NSEG(n)-1)/(NSEG(n)-1), 2*acos(filD(5:13:end)), 'Color', 'k', 'Marker', NSEG_line_symbols{n}, 'DisplayName', sprintf('$\\Delta L/L = %.3f$', D(n,1)));
%     
%     for k=1:NSEG(n)-1
%         
%         theta_sim = 2*acos(filD(5 + 13*(k-1)));
%         
%         left_val = (theta(k) - theta_sim)^2;
%         
%         theta_sim = 2*acos(filD(5 + 13*k));
%         
%         right_val = (theta(k+1) - theta_sim)^2;
%         
%         D(n,2) = D(n,2) + 0.5*(left_val + right_val);
%         
%     end
%     
%     D(n,2) = sqrt(D(n,2)*L/NSEG(n));
%     
% end
% 
% plot((0 : NSEG(n)-1)/(NSEG(n)-1), theta, 'k--', 'DisplayName', 'Exact');
% set(gca,'FontSize',24,'FontName','Times');
% xlabel('$s/L$', 'Interpreter', 'latex');
% ylabel('$\theta(s)$', 'Interpreter', 'latex');
% yticks([4 5 6 7 8]*pi/16);
% yticklabels({'\pi/4', '5\pi/16', '6\pi/16', '7\pi/16', '\pi/2'});
% legend1 = legend(gca,'show');
% legend1.Interpreter = 'latex';
% axis tight;
% 
% % Plot the errors
% 
% figure;
% 
% loglog(D(:,1), D(:,2), 'ko-', 'DisplayName', 'Simulation');
% 
% set(gca,'FontSize',24,'FontName','Times');
% xlabel('$\Delta L/L$','Interpreter','latex');
% ylabel('$E_2$','Interpreter','latex');
% 
% hold on;
% 
% loglog(linspace(D(1,1), D(end,1)), (linspace(D(1,1), D(end,1)).^1), 'k--', 'DisplayName', '$\mathcal{O}(\Delta L)$');
% loglog(linspace(D(1,1), D(end,1)), (linspace(D(1,1), D(end,1)).^2), 'k:', 'DisplayName', '$\mathcal{O}(\Delta L^2)$');
% 
% legend1 = legend(gca,'show');
% legend1.Interpreter = 'latex';
% 























% Finding theta(L)
myfun = @(thetaL) integral(@(x) 1./(sqrt(cos(thetaL)-cos(x))),thetaL,pi/2);

%fun = @(thetaL) sqrt(Bend/(2*F))*myfun(thetaL)-L;
fun = @(thetaL) sqrt(500)*myfun(thetaL)-L;

thetaL = fzero(fun, [0.1 pi/2-0.1])
%thetaL = 8.054788116455077e-01;
%thetaL = 0.858441596177492;
%thetaL=    0.858393711901002
%thetaL = 0.853652402681338;

% Finding theta(s) at each s
theta = zeros(N,1);

theta(N) = thetaL;
theta(1)=pi/2;

i = N-1;

 
for i = N-1 : -1 : 2

    s = (i-1) * L/(N-1);
    
    myfun = @(thetas) integral(@(x) 1./(sqrt(cos(thetaL)-cos(x))),thetas,pi/2);
    
    fun = @(thetas) sqrt(Bend/(2*F))*myfun(thetas)-s;

    theta(i) = fzero(fun, [thetaL pi/2]);

end
theta

% Theta from downwards force code 

theta_code = [1.570796326794896; 1.493780490968601; 1.421574922917614; 
    1.354141997759587; 1.291422434485647; 1.233341677430593; 1.179814350469897;
    1.130748275910069; 1.086047856205970; 1.045616973377673; 1.009361349905147;
    0.977190420369966; 0.949019139022074; 0.924769801240838; 0.904373327521892;
    0.887769802546177; 0.874908939444139; 0.865750981161033; 0.860267521099117;
    0.858441596177492]; %dL = 2.2

theta_code = [1.570796326794896;
   1.493780803822800;
   1.421576486366608;
   1.354144363870832;
   1.291424645620452;
   1.233343276320712;
   1.179815693149450;
   1.130750107065499;
   1.086050672166536;
   1.045620633901879;
   1.009365216127505;
   0.977193953918501;
   0.949022362635081;
   0.924773097990348;
   0.904376903498406;
   0.887773552256024;
   0.874912769027026;
   0.865754972349045;
   0.860271716963047;
   0.858445877473401];

% theta =
% 
%    1.570796326794896
%    1.493779487975762
%    1.421575723817241
%    1.354144224224131
%    1.291425532731375
%    1.233345172493594
%    1.179818077462201
%    1.130752294979963
%    1.086052229973461
%    1.045621572350345
%    1.009365863973955
%    0.977194638892696
%    0.949023153590222
%    0.924773792594500
%    0.904377230453842
%    0.887773388681978
%    0.874912198006625
%    0.865754175415791
%    0.860270837341106
%    0.858444975635821


% Calculating error

error=0;

for i=1:N
   
    error = error + (theta_code(i) - theta(i))^2;

end

error = sqrt(dL*error)


%%
% Plotting together
Ytheory = zeros(3,N);
Ysimulation = zeros(3,N);

for i=2:N
    
    tn = [cos(theta(i));0;sin(theta(i))];
    tn1 = [cos(theta(i-1));0;sin(theta(i-1))];
    
    Ytheory(:,i) = Ytheory(:,i-1) + 0.5*dL*(tn + tn1);

    tns = [cos(theta_code(i));0;sin(theta_code(i))];
    tns1 =  [cos(theta_code(i-1));0;sin(theta_code(i-1))];

    Ysimulation(:,i) = Ysimulation(:,i-1) + 0.5*dL*(tns + tns1);


end

figure;hold on;

plot3([0 Ytheory(1,1)],[0 Ytheory(2,1)],[0 Ytheory(3,1)],'k-o', 'LineWidth',3, 'MarkerSize',10);

plot3([0 Ysimulation(1,1)],[0 Ysimulation(2,1)],[0 Ysimulation(3,1)],'r-o', 'LineWidth',3, 'MarkerSize',10);

for i=2:N
    
    plot3([Ysimulation(1,i-1) Ysimulation(1,i)],[Ysimulation(2,i-1) Ysimulation(2,i)],[Ysimulation(3,i-1) Ysimulation(3,i)],'r-o', 'LineWidth',3, 'MarkerSize',10);

    plot3([Ytheory(1,i-1) Ytheory(1,i)],[Ytheory(2,i-1) Ytheory(2,i)],[Ytheory(3,i-1) Ytheory(3,i)],'k-o', 'LineWidth',3, 'MarkerSize',10);

end

grid on

zlim([-1, 49]); xlim([-30,30]); ylim([-30,30]); 

xlabel('X'); ylabel('Y'); zlabel('Z')

view(0,0)





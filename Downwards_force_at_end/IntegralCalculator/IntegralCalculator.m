% Solving BVP
% ode = (diff(theta,s))^2 == 2*F/Bend * (A - cos(theta));
% cond = theta(0) == pi/2;
% thetaSol(s) = dsolve(ode,cond)


% sspan = [0 20];
% theta0 = pi/2;
% [s,theta] = ode45(@(s,theta) sqrt(2*F/Bend * (cos(5) - cos(theta))), sspan, theta0)



smesh = linspace(1,20,20);
solinit = bvpinit(smesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);
%plot(sol.x, sol.y, '-o')


% Theta from downwards force code
theta_code = [1.570796326794896; 1.493780490968601; 1.421574922917614; 
    1.354141997759587; 1.291422434485647; 1.233341677430593; 1.179814350469897;
    1.130748275910069; 1.086047856205970; 1.045616973377673; 1.009361349905147;
    0.977190420369966; 0.949019139022074; 0.924769801240838; 0.904373327521892;
    0.887769802546177; 0.874908939444139; 0.865750981161033; 0.860267521099117;
    0.858441596177492];

dL = 2.2;

% Calculating error
error=0;
for i=1:length(theta_code)
   
    error = error + (theta_code(i) - (sol.y(1,i)))^2;

end

error = sqrt(dL*error)


% Plotting the filaments to compare
t= sol.y(1,:);
N=20;
Y = zeros(3,N);
Ysim = zeros(3,N);

for i=2:N
    
    tn = [cos(t(i));0;sin(t(i))];
    tn1 = [cos(t(i-1));0;sin(t(i-1))];
    
    Y(:,i) = Y(:,i-1) + 0.5*dL*(tn + tn1);

    tns = [cos(theta_code(i));0;sin(theta_code(i))];
    tns1 =  [cos(theta_code(i-1));0;sin(theta_code(i-1))];

    Ysim(:,i) = Ysim(:,i-1) + 0.5*dL*(tns + tns1);


end

figure;hold on;

plot3([0 Y(1,1)],[0 Y(2,1)],[0 Y(3,1)],'k-o', 'LineWidth',3, 'MarkerSize',10);
plot3([0 Ysim(1,1)],[0 Ysim(2,1)],[0 Ysim(3,1)],'r-o', 'LineWidth',3, 'MarkerSize',10);
for i=2:N
    
    plot3([Y(1,i-1) Y(1,i)],[Y(2,i-1) Y(2,i)],[Y(3,i-1) Y(3,i)],'k-o', 'LineWidth',3, 'MarkerSize',10);
    plot3([Ysim(1,i-1) Ysim(1,i)],[Ysim(2,i-1) Ysim(2,i)],[Ysim(3,i-1) Ysim(3,i)],'r-o', 'LineWidth',3, 'MarkerSize',10);

end



grid on
zlim([-1, 49]);
xlim([-30,30]);
ylim([-30,30]); 

xlabel('X')
ylabel('Y')
zlabel('Z')
view(0,0)

%Y_{n+1} = Y_{n} + 0.5*DL*t_{n} + t_{n+1})


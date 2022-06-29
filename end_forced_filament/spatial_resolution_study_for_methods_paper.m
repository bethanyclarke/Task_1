
% This script compares the centreline orientation of a clamped filament
% bent by a point force applied at the free end to the analytical solution
% given on p. 74 of 'Theory of Elasticity' (3rd edition) by Landau and
% Lifschitz.

clear
close all
format long 

NSEG = [5 10 20 40 80 160];

NSEG_line_symbols = {'o', 's', '^', '<', '>', '.'};

D = zeros(length(NSEG), 2);

L = 44;
f_fac = 1/sqrt(500); % 0.5*KB/f = 0.5*KB/(0.001*KB) = 500.

integral = @(theta0, theta) abs(real((2/sqrt(cos(theta0)-1))*(ellipticF(pi/4, -2/(cos(theta0)-1)) - ellipticF(0.5*theta, -2/(cos(theta0)-1)))));

D(:,1) = 1./NSEG';

figure;
hold on;

for n=1:length(NSEG)
    
    theta_file_name = sprintf('%i_segment_filament_theta_values.dat', NSEG(n));
    
    if exist(theta_file_name, 'file')
        
        theta = load(theta_file_name);
        
        
    else
        
        fprintf('No theta data file was found for a %i-segment filament. Calculating analytical values...\n', NSEG(n));
        
        theta0 = 8.054788116455077e-01; % Previously calculated solution to integral(theta0, theta0) = L/sqrt(0.5*KB/f) = 44/sqrt(500).
        
        theta(NSEG(n)) = theta0; % theta(k) = \theta(s = (k-1)*L/(NSEG(n)-1)).
        theta(1) = pi/2;
        
        for k=2:NSEG(n)-1
            
            upper_bound = theta(k-1);
            lower_bound = theta0;
            
            theta(k) = 0.5*(upper_bound + lower_bound);
            
            s = (k-1)*L/(NSEG(n)-1);
            
            error = s*f_fac - integral(theta0, theta(k));
            
            while abs(error) > 1e-9
                
                if error > 0
                    
                    upper_bound = theta(k);
                    
                else
                    
                    lower_bound = theta(k);
                    
                end
                
                theta(k) = 0.5*(upper_bound + lower_bound);
                
                error = s*f_fac - integral(theta0, theta(k));
                
            end
            
            fprintf('Calculated theta(%i) = %g\n', k, theta(k));
            
        end
        
        dlmwrite(theta_file_name, theta, 'precision', '%.12e');
        
    end
    
    filD = load(sprintf('%isegs.dat', NSEG(n)));
    
    plot((0 : NSEG(n)-1)/(NSEG(n)-1), 2*acos(filD(5:13:end)), 'Color', 'k', 'Marker', NSEG_line_symbols{n}, 'DisplayName', sprintf('$\\Delta L/L = %.3f$', D(n,1)));
    
    for k=1:NSEG(n)-1
        
        theta_sim = 2*acos(filD(5 + 13*(k-1)));
        
        left_val = (theta(k) - theta_sim)^2;
        
        theta_sim = 2*acos(filD(5 + 13*k));
        
        right_val = (theta(k+1) - theta_sim)^2;
        
        D(n,2) = D(n,2) + 0.5*(left_val + right_val);
        
    end
    
    D(n,2) = sqrt(D(n,2)*L/NSEG(n));
    
end

plot((0 : NSEG(n)-1)/(NSEG(n)-1), theta, 'k--', 'DisplayName', 'Exact');
set(gca,'FontSize',24,'FontName','Times');
xlabel('$s/L$', 'Interpreter', 'latex');
ylabel('$\theta(s)$', 'Interpreter', 'latex');
yticks([4 5 6 7 8]*pi/16);
yticklabels({'\pi/4', '5\pi/16', '6\pi/16', '7\pi/16', '\pi/2'});
legend1 = legend(gca,'show');
legend1.Interpreter = 'latex';
axis tight;

% Plot the errors

figure;

loglog(D(:,1), D(:,2), 'ko-', 'DisplayName', 'Simulation');

set(gca,'FontSize',24,'FontName','Times');
xlabel('$\Delta L/L$','Interpreter','latex');
ylabel('$E_2$','Interpreter','latex');

hold on;

loglog(linspace(D(1,1), D(end,1)), (linspace(D(1,1), D(end,1)).^1), 'k--', 'DisplayName', '$\mathcal{O}(\Delta L)$');
loglog(linspace(D(1,1), D(end,1)), (linspace(D(1,1), D(end,1)).^2), 'k:', 'DisplayName', '$\mathcal{O}(\Delta L^2)$');

legend1 = legend(gca,'show');
legend1.Interpreter = 'latex';


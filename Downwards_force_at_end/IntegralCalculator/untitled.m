%% Using Tims code for theory (looks at elliptical integrals, uses bisection method)

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

% Generate theory theta
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

  %  plot((0 : NSEG(n)-1)/(NSEG(n)-1), 2*acos(filD(5:13:end)), 'Color', 'k', 'Marker', NSEG_line_symbols{n}, 'DisplayName', sprintf('$\\Delta L/L = %.3f$', D(n,1)));
   
    % Writing the simulation theta
    
    switch NSEG(n)
        case 5
            theta_code = [1.570796326794896; 1.221309038622039;
                0.985716858705273; 0.851043575734326; 0.807346466359424];
        case 10
            theta_code = [1.570796326794896;
   1.400710730599122;
   1.254197681488358;
   1.130406615491084;
   1.028243269772683;
   0.946552854931199;
   0.884256054407150;
   0.840445786608591;
   0.814441279253425;
   0.805820622124189];
        case 20
            theta_code = [1.570796326794896;
   1.487270447646264;
   1.409090104024004;
   1.336203599210782;
   1.268533562707678;
   1.205983676141971;
   1.148444171872750;
   1.095796713169530;
   1.047918774149257;
   1.004687429804518;
   0.965982534454400;
   0.931689346940301;
   0.901700680098634;
   0.875918644306856;
   0.854256047683248;
   0.836637507848006;
   0.823000318449586;
   0.813295103259450;
   0.807486286351247;
   0.805552406060519];
        case 40
            theta_code = [1.570796326794896
   1.529438281334394
   1.489350424735836
   1.450531014525101
   1.412975400756163
   1.376676880508931
   1.341627527371846
   1.307818064455589
   1.275237818129067
   1.243875408925089
   1.213718329751996
   1.184753604914526
   1.156967926685780
   1.130347406500042
   1.104878283452987
   1.080546468817110
   1.057337770391265
   1.035238173962843
   1.014233415275259
   0.994309568017417
   0.975452902936132
   0.957649885902568
   0.940887673914624
   0.925153771841439
   0.910436307600947
   0.896724108500406
   0.884006421132651
   0.872273194695964
   0.861514927860987
   0.851722666987576
   0.842888206587178
   0.835003984235255
   0.828063221069965
   0.822059969637882
   0.816989037163314
   0.812846086172373
   0.809627540468585
   0.807330598488402
   0.805953254047306
   0.805494277996525];
        case 80

            theta_code = [1.570796326794896;
   1.550219665990420;
   1.529953147075992;
   1.509996579909204;
   1.490349650231679;
   1.471011927931959;
   1.451982872988131;
   1.433261841329495;
   1.414848090485506;
   1.396740785119196;
   1.378939002483568;
   1.361441737756303;
   1.344247909272199;
   1.327356363633085;
   1.310765880703966;
   1.294475178482561;
   1.278482917844296;
   1.262787707158802;
   1.247388106779149;
   1.232282633397203;
   1.217469764270066;
   1.202947941313359;
   1.188715575064930;
   1.174771048515582;
   1.161112720809720;
   1.147738930821081;
   1.134648000593957;
   1.121838238663619;
   1.109307943251476;
   1.097055405338889;
   1.085078911619071;
   1.073376747334376;
   1.061947198997154;
   1.050788556999985;
   1.039899118116586;
   1.029277187890585;
   1.018921082932467;
   1.008829133104312;
   0.998999683613979;
   0.989431097009026;
   0.980121755083799;
   0.971070060694715;
   0.962274439489869;
   0.953733341557109;
   0.945445242989777;
   0.937408647376350;
   0.929622087215957;
   0.922084125261171;
   0.914793355791843;
   0.907748405822922;
   0.900947936248203;
   0.894390642923071;
   0.888075257686497;
   0.882000549328791;
   0.876165324503532;
   0.870568428587819;
   0.865208746492499;
   0.860085203424043;
   0.855196765601253;
   0.850542440926237;
   0.846121279615268;
   0.841932374787133;
   0.837974863013247;
   0.834247924830532;
   0.830750785217890;
   0.827482714038271;
   0.824443026447735;
   0.821631083271674;
   0.819046291351803;
   0.816688103860785;
   0.814556020590295;
   0.812649588210037;
   0.810968400499616;
   0.809512098554463;
   0.808280370965542;
   0.807272953974099;
   0.806489631601951;
   0.805930235756947;
   0.805594646315459;
   0.805482791180324];

        case 160
            theta_code = [1.570796326794896;
   1.560533968663006;
   1.550348186801546;
   1.540238968435971;
   1.530206293845326;
   1.520250135652374;
   1.510370459177050;
   1.500567222485074;
   1.490840376582700;
   1.481189865580524;
   1.471615626886333;
   1.462117591378593;
   1.452695683580556;
   1.443349821834418;
   1.434079918485085;
   1.424885880023513;
   1.415767607274618;
   1.406724995559465;
   1.397757934855561;
   1.388866309961756;
   1.380050000657826;
   1.371308881862926;
   1.362642823795370;
   1.354051692126319;
   1.345535348133995;
   1.337093648854387;
   1.328726447233876;
   1.320433592280153;
   1.312214929202887;
   1.304070299563967;
   1.295999541417346;
   1.288002489452535;
   1.280078975133159;
   1.272228826833705;
   1.264451869975984;
   1.256747927161261;
   1.249116818302896;
   1.241558360755118;
   1.234072369440564;
   1.226658656975392;
   1.219317033794232;
   1.212047308268756;
   1.204849286831641;
   1.197722774088387;
   1.190667572938242;
   1.183683484683850;
   1.176770309144836;
   1.169927844765212;
   1.163155888729445;
   1.156454237054556;
   1.149822684703603;
   1.143261025685510;
   1.136769053154939;
   1.130346559510485;
   1.123993336490815;
   1.117709175269920;
   1.111493866549523;
   1.105347200649704;
   1.099268967598343;
   1.093258957218377;
   1.087316959213716;
   1.081442763253235;
   1.075636159051838;
   1.069896936452319;
   1.064224885503299;
   1.058619796537240;
   1.053081460244833;
   1.047609667750712;
   1.042204210685799;
   1.036864881250607;
   1.031591472297522;
   1.026383777389369;
   1.021241590868236;
   1.016164707919941;
   1.011152924637184;
   1.006206038080871;
   1.001323846343397;
   0.996506148603047;
   0.991752745184833;
   0.987063437616033;
   0.982438028681087;
   0.977876322475315;
   0.973378124456977;
   0.968943241498899;
   0.964571481938079;
   0.960262655624260;
   0.956016573967418;
   0.951833049984494;
   0.947711898344022;
   0.943652935410290;
   0.939655979286274;
   0.935720849855534;
   0.931847368822504;
   0.928035359752791;
   0.924284648111349;
   0.920595061300285;
   0.916966428695685;
   0.913398581682921;
   0.909891353691775;
   0.906444580229842;
   0.903058098915763;
   0.899731749510359;
   0.896465373949360;
   0.893258816371990;
   0.890111923151116;
   0.887024542921853;
   0.883996526608791;
   0.881027727453637;
   0.878118001040839;
   0.875267205323165;
   0.872475200646414;
   0.869741849773287;
   0.867067017906528;
   0.864450572711578;
   0.861892384338580;
   0.859392325442894;
   0.856950271206180;
   0.854566099355831;
   0.852239690184521;
   0.849970926568531;
   0.847759693986044;
   0.845605880533906;
   0.843509376945286;
   0.841470076605118;
   0.839487875566147;
   0.837562672564175;
   0.835694369032233;
   0.833882869115153;
   0.832128079682493;
   0.830429910342106;
   0.828788273452251;
   0.827203084134131;
   0.825674260282791;
   0.824201722578626;
   0.822785394498100;
   0.821425202323618;
   0.820121075153431;
   0.818872944910880;
   0.817680746353054;
   0.816544417079418;
   0.815463897539937;
   0.814439131042250;
   0.813470063759199;
   0.812556644735491;
   0.811698825894118;
   0.810896562042160;
   0.810149810876767;
   0.809458532989851;
   0.808822691873390;
   0.808242253923796;
   0.807717188445761;
   0.807247467656356;
   0.806833066688023;
   0.806473963591833;
   0.806170139340063;
   0.805921577828467;
   0.805728265878319;
   0.805590193237687;
   0.805507352582976;
   0.805479739519737];
          
        otherwise
        
    end

    for k=1:NSEG(n)-1
                
        theta_sim = theta_code(k);

        left_val = (theta(k) - theta_sim)^2;
        
        %theta_sim = 2*acos(filD(5 + 13*k));
        
        %right_val = (theta(k+1) - theta_sim)^2;
        
        D(n,2) = D(n,2) + left_val;
        
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




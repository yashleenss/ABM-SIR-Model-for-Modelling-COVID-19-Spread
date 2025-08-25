%-------------------------------------------------------------------------%
                          % INITIALIZATION %
%-------------------------------------------------------------------------%

  clear

% Number of particles
  N = 5000;                      % total number of Particles
  NRc = 0;                       % initially no one is recovered

% Age Group distribution (percentages)
  elderly_percentage = 0.17;     % 17% are elderly
  adult_percentage = 0.63;       % 63% are adults 
  kids_percentage = 0.20;        % 20% are kids  

% Number of individuals in each group
  n_elderly = round(N * elderly_percentage);
  n_adults = round(N * adult_percentage);
  n_kids = N - n_elderly - n_adults; % remaining are kids

% Initialize age group for tracking the daily infection cases
  age_group = zeros(N, 1); 
  age_group(1:n_elderly) = 1;
  age_group(n_elderly+1:n_elderly+n_adults) = 2;
  age_group(n_elderly+n_adults+1:end) = 3;
  
%---------------------------- SIMULATION AREA ----------------------------%

% Initialize the length of the simulation area
  density = 0.663;          % initial density (particles/m^2)
  area = N / density;
  L = sqrt(area);

% Initialize particle positions (random initial positions within the Lo)
  X = L * rand(N, 1);
  Y = L * rand(N, 1);

%--------------------------- TIME SPAN VARIABLES -------------------------% 

% Number of simulation days (for example first outbreak = 90 days)
  T_span = 500;              

% Distinct theoretical reproduction number for each outbreak
  R_t_vector = zeros(T_span,1);
  R_t_vector(1:199) = 1.168293844;   % First outbreak (decay phase)
  R_t_vector(200:364) = 2.497643506; % Second outbreak (delta variant)
  R_t_vector(365:end) = 2.5;         % Third outbreak (omicron variant)  

% Preallocate daily infected cases per age group
  daily_new_kids = zeros(1, T_span); 
  daily_new_adults = zeros(1, T_span);
  daily_new_elderly = zeros(1, T_span);

% Preallocate the sum of daily infected cases across simulations
  sum_daily_kids = zeros(1, T_span);
  sum_daily_adults = zeros(1, T_span);
  sum_daily_elderly = zeros(1, T_span);

% Preallocate the sum of average infection probability across simulations
  sum_avg_kids = zeros(1, T_span);
  sum_avg_adults = zeros(1, T_span);
  sum_avg_elderly = zeros(1, T_span);
  
% Preallocate the sum of health status and infection probability status
  sum_H = zeros(N,1);
  sum_P_H = zeros(N, 1);

% Preallocate the sum of no. of individuals in each state
  sum_infected = zeros(1,T_span);
  sum_recovered = zeros(1,T_span);
  sum_susceptible = zeros(1, T_span);

% Preallocate to store the VE vectors in each state and variant
  sum_VE_all = zeros(1, T_span);
  sum_VE_infected = zeros(1, T_span);
  sum_VE_recovery = zeros(1, T_span);
  sum_VE_susceptible = zeros(1, T_span);
  
%------------------------- VACCINE DISTRIBUTION --------------------------%

% Initialize the change in vaccination through VEdays and  VE levels
  VElevel = 0.9;                  % maximum VElevel is 90%
  VEdays = [7, 28, 150, 180];     % days after which VE changes

% Vaccination coverage for each age group
  coverage_kids = 0.60;
  coverage_adults = 0.80;
  coverage_elderly = 0.95;

% Number of people to be vaccinated in each age group
  num_elderly_vaccinated = round(n_elderly * coverage_elderly);
  num_adults_vaccinated = round(n_adults * coverage_adults);
  num_kids_vaccinated = round(n_kids * coverage_kids);

% Unique recovery time for each individual
  mu1 = 14; sigma1 = 4;
  mu2 = 18; sigma2 = 6;

% Recovery time averaging around 16, between [7, 21] days
  raw = (normrnd(mu1, sigma1, [N,1])/2 + normrnd(mu2, sigma2, [N,1])/2);
  RecoveryTime = ceil(min(max(raw, 7), 21));

% Adjusting immunity duration wrt. time span, bounding to [90, 180] days 
  ImmunityTime = ceil(min(max(raw + 100, 90), 180));

% Preallocating variables for checking vaccination stratified behavior
  ngps = 3;                         % define VE based number of groups
  VE_fixed = [0.3, 0.7];            % vaccinantion cut-offs
  VE_group = zeros(N,1);             
  P_H_tspan = zeros(ngps, T_span);
  sum_P_H_tspan = zeros(ngps, T_span);

%------------------------------ RADII VARIABLES --------------------------%

% Infection radius to control individual's mobility
  infection_radii = 6;               
  radii_squared = infection_radii.^2;
  num_radii = length(infection_radii); % Number of radii

% Initialize storage for mean infection probability over time per radius
  Infected_over_radii = zeros(T_span, num_radii);
  Mean_PH_over_radii = zeros(T_span, num_radii); 

% Preallocate the array to avoid resizing
  plot_handles = gobjects(num_radii, 1);  

%-------------------------------- VIDEO WRITER ---------------------------%

% Create an object for animation
  videoFilename = VideoWriter('particles_animation.mp4', 'MPEG-4');
  v = (videoFilename);
  v.FrameRate = 60;          % Set frame rate for better resolution
  open(v);

%----------------- RADIUS LOOP -------------------%

% FOR Loop over radii (0:6)
  for r_idx = 1:num_radii
      radius = radii_squared(r_idx); % Current radius
  
%--------------- SIMULATION LOOP ------------------%  

% Define the number of simulations
  num_simulations = 500;
  for sim = 1:num_simulations

      %--------------------------- HEALTH STATUS -------------------------%

      % Re-initialize health status: 1 for infected, 0 for healthy
        [H] = HealthStatus(N, NRc);
        infected_idx = find(H == 1);
        susceptible_idx = find(H == 0);
        
      % Re-assign vaccination days for each age group
        Tv = Inf(N, 3);   % initialize Tv

      % First vaccination campaign Day 28 to Day 70
        vaccinated_elderly_1 = randsample(find(age_group == 1),...
                                               num_elderly_vaccinated);
        Tv(vaccinated_elderly_1,1) = randi([28, 48],...
                                               num_elderly_vaccinated, 1);
        vaccinated_adults_1 = randsample(find(age_group == 2),...
                                                num_adults_vaccinated);
        Tv(vaccinated_adults_1,1) = randi([41, 70],...
                                                num_adults_vaccinated, 1);
        Tv(n_elderly+n_adults+1:end, 1) = Inf; % kids are not vaccinated

      % Second vaccination campaign Day 91 to Day 207
        vaccinated_elderly_2 = randsample(find(age_group == 1),...
                                               num_elderly_vaccinated);
        Tv(vaccinated_elderly_2, 2) = randi([91, 112],...
                                               num_elderly_vaccinated, 1);
        vaccinated_adults_2 = randsample(find(age_group == 2),...
                                               num_adults_vaccinated);
        Tv(vaccinated_adults_2, 2) = randi([113, 149],...
                                               num_adults_vaccinated, 1);
        vaccinated_kids_2 = randsample(find(age_group == 3),...
                                               num_kids_vaccinated);
        Tv(vaccinated_kids_2, 2) = randi([150, 207],...
                                               num_kids_vaccinated, 1);
      
      % Third vaccination campaign Day 232 to Day 303
        vaccinated_elderly_3 = randsample(find(age_group == 1),...
                                               num_elderly_vaccinated);
        Tv(vaccinated_elderly_3, 3) = randi([232, 253],...
                                               num_elderly_vaccinated, 1);
        vaccinated_adults_3 = randsample(find(age_group == 2),...
                                               num_adults_vaccinated);
        Tv(vaccinated_adults_3, 3) = randi([254, 284],...
                                               num_adults_vaccinated, 1);
        vaccinated_kids_3 = randsample(find(age_group == 3),...
                                               num_kids_vaccinated);
        Tv(vaccinated_kids_3, 3) = randi([285, 303],...
                                               num_kids_vaccinated, 1);

        Ti = Tv;                          % initialized recovery 
        VE = zeros(N, 1);                 % reinitialize at each run
        Immunity = zeros(N,1);            % vaccine = 0 & recovery = 1

      %----------------------- INFECTION STATE ---------------------------%

      % Initialize tI & random infection days ONLY to infected individuals
        tI = zeros(N,1);
        infection = zeros(1,N);
        tI(infected_idx) = randi([0,500], length(infected_idx),1);

      %--------------------- EXTRACTING RESULTS --------------------------%
      
       results = MAM(N, L, H, X, Y, T_span, Tv, Ti, tI, VE, VElevel,...
                     VEdays, VE_fixed, ngps, age_group, infection,...
                     radius, RecoveryTime, Immunity, ImmunityTime,...
                     daily_new_elderly, R_t_vector, P_H_tspan,...
                     daily_new_adults, daily_new_kids);

      % 4.2 - Accumulate the daily infected cases for S,I,R dynamics
        sum_infected = sum_infected + results.daily_infected_idx;
        sum_recovered = sum_recovered + results.daily_recovery_idx;
        sum_susceptible = sum_susceptible + results.daily_susceptible_idx;
 
      % 4.3.1 - Accumulate the total daily infections per day per age group
        sum_daily_kids = sum_daily_kids + results.daily_new_kids;
        sum_daily_adults = sum_daily_adults + results.daily_new_adults;
        sum_daily_elderly = sum_daily_elderly + results.daily_new_elderly;
  
      % 4.3.2 - Accumulate the average infection probability per age group
        sum_avg_kids = sum_avg_kids + results.Avg.kids;
        sum_avg_adults = sum_avg_adults + results.Avg.adults;
        sum_avg_elderly = sum_avg_elderly + results.Avg.elderly;
  
      % 4.4 - Accumulate the total health status and infection probability
        sum_H = sum_H + results.H;
        sum_P_H = sum_P_H + results.P_H;
        sum_P_H_tspan = sum_P_H_tspan + results.P_H_tspan;
 
      % 4.5 - Accumulate the S,I,R dynamics for vaccine effectiveness
        sum_VE_infected = sum_VE_infected + results.VE_infected;
        sum_VE_recovery = sum_VE_recovery + results.VE_recovery;
        sum_VE_susceptible = sum_VE_susceptible + results.VE_susceptible;
 
      % 4.6 -Accumulate the variant type dynamics for vaccine effectiveness
        sum_VE_all = sum_VE_all + results.VE_all;
      
  end % end simulation loop

% 4.7 - Update results for infected individuals for the present radius
  Infected_over_radii = sum_infected(:) / num_simulations;

% Normalize by total population to get mean infection probability 
  Mean_PH_over_radii(:, r_idx) = Infected_over_radii / N;

  end % end radii loop

%-------------------------------------------------------------------------%
%                             ANALYSIS TESTS
%-------------------------------------------------------------------------%

% [R G B] Color Palette 
  light_blue = [0.1, 0.44, 0.9];
  light_red = [0.95, 0.07, 0.09];
  light_green = [0.37, 0.9, 0.18];

% Test 4.2: Total number of susceptible, infected & recovered individuals
  figure;
  plot(1:T_span, sum_susceptible/num_simulations, '-', 'Color',...
                 light_blue, 'LineWidth', 2, 'MarkerSize', 3,...
                 'MarkerFaceColor', light_blue); % Susceptible
  hold on;
  plot(1:T_span, sum_infected/num_simulations, '-', 'Color',...
                 light_red, 'LineWidth', 2, 'MarkerSize', 3,...
                 'MarkerFaceColor', light_red); % Infected

  plot(1:T_span, sum_recovered/num_simulations, '-', 'Color',...
                 light_green, 'LineWidth', 2, 'MarkerSize', 5,...
                 'MarkerFaceColor', light_green); % Infected

% Highlight the infection peak
  [max_infected, max_time] = max(sum_infected / num_simulations);
  plot(max_time, max_infected, 'p', 'MarkerSize', 10, ...
                 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow');

% Add the infection peak label
  text(max_time + 2, max_infected,sprintf('Peak: %d',...
                     round(max_infected)),'FontSize', 10,...
                     'FontWeight', 'bold','Color', 'black');

 % xlabel('Time (days)','FontSize', 12);
  ylabel('Mean Number of Individuals per state (Radius = 6 cm)',...
                                                'FontSize', 12);
  title(['Temporal Dynamics of COVID-19 Transmission across 500 ' ...
                                 'simulations'], 'FontSize', 12);
  legend({'Susceptible', 'Infected','Recovered'},'Location',...
                                'northeast', 'FontSize', 10);
  set(gca, 'FontSize', 10, 'LineWidth', 1);

% Custom X-axis design according to date ticks
  start_date = datetime(2020,12,12);
  date_ticks = start_date + days([0 27 90 231 500]);
  xticks([1 28 91 232 500])
  xticklabels(datestr(date_ticks, 'mmm dd, yyyy'))

% Add vertical lines on major dates
  xline(28, '--k','LineWidth', 2,'HandleVisibility', 'off');
  xline(91, '--k','LineWidth', 2,'HandleVisibility','off');
  xline(232, '--k','LineWidth', 2,'HandleVisibility','off');
  yLimits = ylim;

% Light Grey Patch for First vaccination campaign
  fill([28 71 71 28], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)],...
                [0.9 0.9 0.9],'FaceAlpha', 0.3, 'EdgeColor', 'none',...
                                          'HandleVisibility','off');
  text(28 + 2, yLimits(2)*0.95, 'First Vaccination Campaign',...
                         'FontSize', 10, 'Color', [0.3 0.3 0.3])
  grid on;
  hold off;

% Test 4.3.1: Daily new infected cases w.r.t. age groups
  figure;
  plot(1:T_span, sum_daily_kids/num_simulations, '-^', 'Color',...
                    light_blue, 'LineWidth', 2, 'MarkerSize', 3,...
                                 'MarkerFaceColor', light_blue); 
  hold on;
  plot(1:T_span, sum_daily_adults/num_simulations, '-o', 'Color',...
                      light_red, 'LineWidth', 2, 'MarkerSize', 3,...
                                   'MarkerFaceColor', light_red); 
  plot(1:T_span, sum_daily_elderly/num_simulations, '-.', 'Color',...
                   light_green, 'LineWidth', 2.5, 'MarkerSize', 5,...
                                  'MarkerFaceColor', light_green); 

  xlabel('Time (days)','FontSize', 12);
  ylabel('Number of New Infections per day per age group (500 runs)',...
                                                     'FontSize', 12);
  title('COVID-19 Daily Confirmed Cases per age group: (2020–2025))',...
                                                     'FontSize', 12);
  legend({'Kids(0-11)', 'Adults (12-59)','Elderly(60+)'},'Location',...
                                        'northeast', 'FontSize', 10);
  set(gca, 'FontSize', 10, 'LineWidth', 1);

  grid on;
  hold off;

% Test 4.3.2: Average P_H per group over Tspan across simulations
  figure;
  plot(1:T_span, sum_avg_elderly/num_simulations, '-o', 'Color',...
                     light_red, 'LineWidth', 2, 'MarkerSize', 3,...
                                  'MarkerFaceColor', light_red);    
  hold on;
  plot(1:T_span, sum_avg_adults/num_simulations, '-', 'Color',...
                    light_green, 'LineWidth', 2, 'MarkerSize', 3,...
                              'MarkerFaceColor', light_green); 
  plot(1:T_span, sum_avg_kids/num_simulations, '-.', 'Color',...
                 light_blue, 'LineWidth', 2, 'MarkerSize', 3,...
                              'MarkerFaceColor', light_blue);  
    
  title('Mean Infection Behavior with three age groups (Radius = 6cm)',...
                                                       'FontSize', 12);
  xlabel('Time (Days)', 'FontSize', 12);
  ylabel('Mean Probability of getting infection (500 simulations)',...
                                                      'FontSize', 12);
  legend({'Elderly (60+)', 'Adults (12-59)','Kids (0-11)'},...
                          'Location', 'best', 'FontSize', 10);     
  set(gca, 'FontSize', 10, 'LineWidth', 1);
  grid on;
  hold off;

% Test 4.4: Mean infection probability by VE stratified groups
  avg_P_H_tspan = sum_P_H_tspan / num_simulations;
  figure; 
  hold on;
  colors = [light_blue; light_red; light_green];

  for i = 1:ngps
      plot(1:T_span, avg_P_H_tspan(i,:), 'Color', colors(i,:),...
                                         'LineWidth', 2);
  end
  xlabel('Time (days)');
  ylabel('Mean Infection Probability of infected population (500 runs)');
  legend({'Group 1 (< 0.3)','Group 2 (0.3-0.7)','Group 3 (≥ 0.7)'});
  title('Mean Infection Behavior with Three VE groups (Radius = 6cm)');
  grid on;
  hold off;

% Test 4.5 - Average VE behavior changes over time per state S,I,R
  Avg_VE_infected = sum_VE_infected / num_simulations;
  Avg_VE_recovery = sum_VE_recovery / num_simulations;
  Avg_VE_susceptible = sum_VE_susceptible / num_simulations;

  figure;
  plot(1:T_span, Avg_VE_susceptible,'-','Color',light_blue,'LineWidth', 2); 
  hold on;
  plot(1:T_span, Avg_VE_infected,'-','Color',light_red,'LineWidth', 2);
  plot(1:T_span, Avg_VE_recovery,'--','Color',light_green,'LineWidth', 2);
  xlabel('Time (days)');
  ylabel('Average Vaccine Effectiveness (radius = 6cm)');
  title('Evoultion of Vaccine Effectiveness Over Time (500 runs)');
  legend({'Susceptible','Infected','Recovered'}, 'Location', 'best');

% Add vertical lines on major dates (3 vaccination campaigns)
  xline(28, '--k','HandleVisibility', 'off');
  xline(91, '--k','HandleVisibility', 'off');
  xline(232, '--k','HandleVisibility', 'off');
  grid on;
  hold off;

% Test 4.6 - Average VE behavior changes over time per variant
% Avg_VE_alphadelta = sum_VE_all / num_simulations; % 60% for Alpha+Delta
  Avg_VE_omicron = sum_VE_all / num_simulations;    % 30% for Omicron

  figure;
  plot(1:T_span, squeeze(Avg_VE_omicron), 'b-', 'LineWidth', 2);
  hold on;
  xlabel('Time (days)');

% Add custom labels on major dates when Time span = 2000 days
  xline(500, '--k','LineWidth', 2);   
  xline(1095, '--r','LineWidth', 2);
  ylabel('Average Vaccine Effectiveness with radius 6cm');
  title('Evoultion of Vaccine Effectiveness Over Time (500 runs)');
  legend({'Omicron Variant','Variants Cutoff',...
          'Post-Pandemic Begins'}, 'Location', 'best');
  grid on;
  hold off;

% Test 4.7 - Number of infected individuals vs Time Steps (for each radii)
  figure;
  hold on;
  colors = lines(num_radii); % Generate distinct colors for each line
  
  for r_idx = 1:num_radii
      
     % Plot infection data for the current radius
       plot_handles(r_idx) = plot(1:T_span,...
                                  Mean_PH_over_radii(:, r_idx),...
                                  '-','Color',...
                                  colors(r_idx, :),'DisplayName',...
                                  sprintf('Radius = %d',...
                                  infection_radii(r_idx)));
  end

  xlabel('Time (days)');
  ylabel('Mean infection probability (500 simulations)');
  title('Infection Probability Over Time for Different Radii (0,1,3,6)cm');
  legend(plot_handles, 'Location', 'Best');
  set(gca,'FontSize',10,'FontName','Times');
  grid on;
  hold off;

% Test 4.8 - Capture the video frame
  frame = getframe(gcf);
  writeVideo(v, frame);
  close(v);            % Close the video file

% Test 4.9 - Final positions and infection status for each particle
  figure;
  gscatter(X, Y, H,[dark_green; dark_red], 'ox');
  title('Individual Agent Position and Infection Status',...
       'FontSize', 12);
  xlabel('X Coordinate');
  ylabel('Y Coordinate');
  axis([0 L 0 L]);
  grid on;
  legend('Susceptible', 'Infected','Location','northeast','FontSize', 12);  

% MAM MODEL
% Initialize variables

% Initialize the length of the simulation area
  Lo = 124;  

% Periodic Boundary length(L) of the simulation area
  L = Lo * (1.25/2.05) ^ (1/2);  % Rt Value for Alpha and Delta Outbreak

% Number of particles
  N = 5000;                       % Total number of Particles
  NRc = 0;                       % Number of intially recovered particles 
% NI = randi([1, N - NRc], 1);   % Total number of infected particles
  NI = round(N * 0.01);          % Starting with 1% infected particles

% Number of days to simulate (for first outbreak = Day(1-90))
  T_span = 90;                   % From Dec.12,2020 to Mar.12,2021  

% Age Group distribution (percentages)
  elderly_percentage = 0.17;     % 17% are elderly
  adult_percentage = 0.63;       % 63% are adults
  kids_percentage = 0.20;        % 20% are kids

% Number of individuals in each group
  n_elderly = round(N * elderly_percentage);
  n_adults = round(N * adult_percentage);
  n_kids = N - n_elderly - n_adults; % Remaining are kids

% Assign age group labels for tracking the daily infection
  age_group = zeros(N, 1);  
  age_group(1:n_elderly) = 1;                        % Elderly first
  age_group(n_elderly+1 : n_elderly+n_adults) = 2;   % Adults next
  age_group(n_elderly+n_adults+1 : end) = 3;         % Kids 

% Preallocate daily cases per age group
  daily_new_elderly = zeros(1, T_span);
  daily_new_adults = zeros(1, T_span);
  daily_new_kids = zeros(1, T_span);

% Standard Deviation
  sigma = Lo / 2; 

% Standard deviation for the infection probability calculation
  sigma_r = sigma/1;

% Initialize health status: 1 for infected, 0 for healthy
  [~, ~, H] = HealthStatus(N, NI, NRc, Lo);
  susceptible_idx = find(H==0);
  infected_idx = find(H==1);

% Select a variant type to work with
  variant = 'Alpha and Delta'; 
  fprintf('Working with %s variant\n', variant);

% Initialize Tv & random vaccination days ONLY to susceptible individuals 
  Tv = zeros(N, 1);  
  Tv(susceptible_idx) = randi([33, 70], length(susceptible_idx),1);

% Assign vaccination days for each group based on first dose (Jan8 - Feb20)
  Tv(1:n_elderly, 1) = randi([1, 10], n_elderly, 1);
  Tv(n_elderly+1:n_elderly+n_adults, 1) = Tv(n_elderly, 1) + 14 + randi([0, 7], n_adults, 1);
  Tv(n_elderly+n_adults+1:end, 1) = 0;
  VE = zeros(N, 1);

% Initialize tI & random infection days ONLY to infected individuals
  tI = zeros(N,1);
  tI(infected_idx) = randi([1, 11], length(infected_idx),1);
  infection = zeros(1,N);

% Estimated infection radius within which particles can infect each other
  infection_radii = 0:6; % inf_rad2 = infection_radii^2;
  num_radii = length(infection_radii); % Number of radii

% Initialize storage for number of infected over time
  num_infected_over_radii = zeros(T_span, num_radii);

% Initialize particle positions (random initial positions within the Lo)
  X = L * rand(N, 1);
  Y = L * rand(N, 1);

%{
% Create a VideoWriter object
  videoFilename = VideoWriter('particles_animation.mp4', 'MPEG-4');
  v = (videoFilename);
  v.FrameRate = 60; % Set frame rate for better resolution
  open(v);
%}

% Initialize storage vectors for time-evolving probabilities 
  P_H0 = zeros(N, T_span);
  active_infections = zeros(1,T_span);

% Preallocate the average infection probability for each radius
  P_H_over_radii = zeros(num_radii, 1); 

% Preallocate the array to avoid resizing
  plot_handles = gobjects(num_radii, 1);

% FOR Loop over radii (0:6)
  for r_idx = 1:num_radii
      radius = infection_radii(r_idx); % Current radius

% Initialize arrays to store the no. of susceptible and infected over time
  num_susceptible = zeros(T_span, 1);
  newly_infected_idx = zeros(T_span, 1);

% Initialize average infection probability over the number of radii 
  sum_prob_sim = zeros(1, T_span);
  
% Defining the total number of simulations
  num_simulations = 1;
  for sim = 1:num_simulations

% Loop through each day to update positions and check infections
  for t = 1:T_span
      
      % Store temporary probabilities for the current time step
        P_H = zeros(N, 1); 

      % Compute daily displacements
        deltaX = normrnd(0, sigma, [N, 1]);
        deltaY = normrnd(0, sigma, [N, 1]);

      % Update New (X,Y) particle positions
        X = mod(X + deltaX, L);
        Y = mod(Y + deltaY, L);

      % Apply periodic boundary conditions (no particle is leaving Lo)
        X(X > L) = X(X > L) - L;
        Y(Y > L) = Y(Y > L) - L; 

      % Check the infection probability for N particles
        for i = 1:N
        
      % If the particle is healthy (H=0)
        if H(i) == 0 
           P_H(i) = 0;   % Define the infection probability
            
           % Check the infection probability for the next particle
             for j = 1:N

             % If the other particle is infected (H=1)
               if i ~= j && H(j) == 1 && (t-tI(j)) > 0 && (t - tI(j)) <= 14

               % Calculate the distance between i and j particles
                 distance = ((X(j) - X(i))^2 + (Y(j) - Y(i))^2); 
                 exp_term = exp(-distance / (2 * (sigma_r)^2));

               % Check the infection radius 
                 if distance <= radius
                    % [t distance inf_rad2]

                    VE(j) = VaccineEffectiveness(Tv(j), t);
                    infection(j) = ComputeInfectious(tI(j), t, variant);
       
                    P_H(i) = P_H(i) + (infection(j)...
                             * (1 - VE(j)) * (exp_term));

                 end 
               end
              end

      % Ensure Clamping values within [0, 1]
        P_H(i) = max(0, min(1, P_H(i))); 

      % Uniformly Distributed number within [0,1] for each particle 
        xi = rand; 

      % Convert to binary (0 or 1)
        if P_H(i) > 0
           P_H(i) = floor(P_H(i) + xi);
        end
                
      % Only update for the NEWLY infected individuals
        if P_H(i) > 0 && H(i) == 0  
           H(i) = 1;     % Mark these particle as infected, H = 1              
           tI(i) = t;    % Set the day of infection for infected particles
                         % No need to check other infected particles
        end
        end
        end
     

% Update the daily number of infected individuals (age groups) over time
  newly_infected_idx = (H == 1) & (tI == t) & (t-tI) <= 14; % (N*1) array vector
  daily_new_elderly(t) = sum(newly_infected_idx & (age_group == 1));
  daily_new_adults(t) = sum(newly_infected_idx & (age_group == 2));
  daily_new_kids(t) = sum(newly_infected_idx & (age_group == 3));

% Calling the avergae function
 % P_H0(:,t)=P_H;
 % Avg = Average(P_H0);

% Total population ever infected until day t
% cumulative_infected(t) = sum(tI > 0 & tI <= t);

% Summation of contagious infected individuals excluding tI=0 noise
  active_infections(t) = sum(H == 1 & tI <= t & (t - tI) <= 14);

% Accumulate mean P_H at each time step t
%  sum_prob_sim(t) = sum_prob_sim(t) + mean(P_H);

  end % end for loop 

  end % end sim loop

% Update results for infected individuals for the present radius
%  num_infected_over_radii(:, r_idx) = num_infected_over_radii(:, r_idx) + num_infected;

% Average across the number of simulations
%  P_H_over_radii(r_idx) = mean(P_H) / sim;
%  num_infected_over_radii(:, r_idx) = num_infected_over_radii(:, r_idx) / sim;

% Average probability over the number of simulations
%  Average_prob_sim = sum_prob_sim / sim;

  end % end radii loop

% RGB Color Palette 
  dark_red = [0.8, 0.1, 0.1];
  dark_green = [0.2, 0.5, 0.2];
  dark_blue = [0.1, 0.3, 0.8];
  light_red = [1, 0.6, 0.6];  

% Test 1: total active infections at every time step
  figure;
 % plot(1:T_span, N * exp(-0.1*(1:T_span)), 'b--');
  plot(1:T_span, max(N) * normpdf(1:T_span, T_span/2, sigma), 'b--');
  hold on;
  plot(1:T_span, active_infections, 'r-', 'LineWidth', 2);
  legend('Normal Distribution','Active Infections','Location','best');
  xlabel('Time (Date)');
  ylabel('Evolution of Active Infections');
  title('Comparison of Normal(Gaussian) Distribution vs Active Infections');

% Add vertical lines on major dates
  xline(28, '--k','HandleVisibility', 'off');
  xline(71, '--k','HandleVisibility', 'off');

% Set custom X-axis ticks
  xticks([1 28 71 90])
  xticklabels({'Dec 20 2020','Jan 8', 'Feb 20','Mar 12,2021'})
  grid on;
  hold off;
% xticks(1:round(T_span/3):T_span)

% Test 2: plot the daily new cases w.r.t. age groups
  figure;
  plot(1:T_span, daily_new_elderly, '-o', 'Color', dark_green, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_green); 
  hold on;
  plot(1:T_span, daily_new_adults, '-*', 'Color', dark_red, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_red); 
  plot(1:T_span, daily_new_kids, '-', 'Color', dark_blue, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_red); 

  xlabel('Time (Date)','FontSize', 12);
  ylabel('Daily New Infected Cases at each time step','FontSize', 12);
  title('Daily Confirmed Cases for First Outbreak ','FontSize', 12);
  legend({'Elderly (60+)', 'Adults (12-59)','Kids (0-11)'},'Location', 'northeast', 'FontSize', 10);
  set(gca, 'FontSize', 10, 'LineWidth', 1);

  % Add vertical lines on major dates
  xline(28, '--k','HandleVisibility', 'off');
  xline(71, '--k','HandleVisibility', 'off');

% Set custom X-axis ticks
  xticks([1 28 71 90])
  xticklabels({'Dec 20 2020','Jan 8', 'Feb 20','Mar 12,2021'})
  grid on;
  hold off;
%{
% Results for Average infection probability w.r.t. radii
  figure;
  plot(infection_radii, P_H_over_radii, '-o', 'LineWidth', 2);
  set(gca,'FontSize',24,'FontName','Times');
  xlabel('Radii');
  ylabel('Average Infection Probability');
  title('Changes in Infection Probability for Different Radii (0 to 6)');
  grid on;

% Results for Infected individuals vs Time Steps (for each radii)
  figure;
  hold on;
  colors = lines(num_radii); % Generate distinct colors for each line
  
  for r_idx = 1:num_radii
      % Plot infection data for the current radius
        plot_handles(r_idx) = plot(1:T_span, num_infected_over_radii(:, r_idx),...
                            'o', 'Color', colors(r_idx, :),...
                            'DisplayName', sprintf('Radius = %d', infection_radii(r_idx)));
  end

  xlabel('Time Span (days)');
  ylabel('Number of Infected Individuals');
  title('Infection Spread Over Time for Different Radii (0 to 6)');
  legend(plot_handles, 'Location', 'Best');
  set(gca,'FontSize',24,'FontName','Times');
  grid on;
  hold off;

% Results for Average Probability over simulation (sim,t)
  figure;
  plot(1:T_span, Average_prob_sim, 'LineWidth', 2);
  set(gca,'FontSize',11,'FontName','Times');
  xlabel('Time (days)');
  ylabel('Average Infection Probability');
  title('Average Infection Probability Over Time');
  grid on; 

% Separate infection probability by Groups
  P_H_group1 = P_H(group1);     
  P_H_group2 = P_H(group2);
  P_H_group3 = P_H(group3);

% Calculate infection proportions for each group
  infected_group1 = mean(P_H_group1); 
  infected_group2 = mean(P_H_group2); 
  infected_group3 = mean(P_H_group3);

% Display infection proportions for each group
  disp(['Proportion Infected (VE >= 90%): ', num2str(infected_group1)]);
  disp(['Proportion Infected (10% <= VE < 90%): ', num2str(infected_group2)]);
  disp(['Proportion Infected (VE < 10%): ', num2str(infected_group3)]);

% Plot Comparison of different VE levels
  figure;
  bar([1, 2, 3], [infected_group1, infected_group2, infected_group3], 'FaceColor', 'flat');
  xticks([1, 2, 3]);
  xticklabels({'VE >= 90%', '10% <= VE < 90%', 'VE < 10%'});
  ylabel('Proportion Infected','FontSize',12);
  title('Comparison of Infection Proportions by Vaccine Effectiveness','FontSize',12);
  grid on;

% RGB Color Palette 
  dark_red = [0.8, 0.1, 0.1];
  dark_green = [0.2, 0.5, 0.2];
  dark_blue = [0.1, 0.3, 0.8];
  light_red = [1, 0.6, 0.6];

% Plot the susceptible and infected individuals over time
  hold on;
  plot(1:T_span, num_susceptible, '-o', 'Color', dark_green, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_green); % Susceptible
  plot(1:T_span, num_infected, '-x', 'Color', dark_red, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_red); % Infected

  xlabel('Time (days)','FontSize', 12);
  ylabel('Individuals at each time step','FontSize', 12);
  title('Evolution of Population Dynamics Over Time','FontSize', 12);
  legend({'Susceptible', 'Infected'},'Location', 'northeast', 'FontSize', 10);
  grid on;
  set(gca, 'FontSize', 10, 'LineWidth', 1);

% Highlight peak of infection 
  [max_infected, max_time] = max(num_infected);
  plot(max_time, max_infected, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow'); % Highlight peak
  text(max_time + 2, max_infected, sprintf('Peak: %d', max_infected), 'FontSize', 10);
  hold off;

% Results of the Average infection/group over time (T_span)
  figure;
  plot(1:T_span, Avg.kids, '-o', 'Color', dark_red, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_red);    % Kids (0-11)
  hold on;
  plot(1:T_span, Avg.adults, '-s', 'Color', dark_green, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_green); % Adults (12-59)
  plot(1:T_span, Avg.elderly, '-^', 'Color', dark_blue, 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', dark_blue);  % Elderly (60+)
    
% Titles and labels
  title('Average Infection Probability Over Time', 'FontSize', 12);
  xlabel('Days', 'FontSize', 12);
  ylabel('Average Infection Probability', 'FontSize', 12);
  legend({'Kids (0-11)', 'Adults (12-59)', 'Elderly (60+)'}, 'Location', 'best', 'FontSize', 10);
  grid on;      
  set(gca, 'FontSize', 10, 'LineWidth', 1);
  hold off;

    
% Capture the frame
  frame = getframe(gcf);
  writeVideo(v, frame);

  end
  
% Close the video file
  close(v);


% Scatter Plot of the final positions and infection status
  figure;
  gscatter(X, Y, H,[dark_green; light_red], 'ox');
  title('Individual Agent Position and Infection Status (OMICRON)','FontSize', 12);
  xlabel('X Coordinate');
  ylabel('Y Coordinate');
  axis([0 L 0 L]);
  grid on;
  legend('Susceptible', 'Infected','Location','Best','FontSize', 10);  
%}
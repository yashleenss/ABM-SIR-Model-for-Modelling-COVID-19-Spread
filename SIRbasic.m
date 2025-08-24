%-------------------------------------------------------------------------%
%           TRADITIONAL Susceptible-Infected-Recovered MODEL
%-------------------------------------------------------------------------%

% Initialize Parameters
  infection_rate = 0.3;               % initial infection rate
  recovery_rate = 0.1;                % initial recovery rate
  R_o = infection_rate / recovery_rate; % basic reproduction number
  T_span = [0 90];                    % simulation duration time            

% Initial conditions
  S0 = 20000;                         % initial susceptible population
  I0 = 1948;                          % initial infected population
  R0 = 0;                             % initial recovered population
  N = S0 + I0 + R0;                   % total population
  initial_conditions = [S0; I0; R0];  % initial conditions for vector
  
% SIR model equations
  SIR_model = @(t, SIR) [
             -infection_rate * SIR(1) * SIR(2) / N;  
             infection_rate * SIR(1) * SIR(2) / N - recovery_rate * SIR(2);
             recovery_rate * SIR(2)                 
                        ];

% Solve the SIR model using ode45
  [t, SIR] = ode45(SIR_model, T_span, initial_conditions);

%---------------------------- RESULT ANALYSIS ----------------------------%
% Subplots of number of individuals in each compartment
  S = SIR(:,1);
  I = SIR(:,2);
  R = SIR(:,3);

% Susceptible - light blue
  subplot(3,1,1);
  area(t, S, 'FaceColor', [0.6 0.8 1], 'EdgeColor', 'none');
  xlabel('Time (days)');
  ylabel('no. of susceptible');
  title('SIR Model - Susceptible Group');
  grid on;

% Infected - Light Red
  subplot(3,1,2);
  area(t, I, 'FaceColor', [1 0.6 0.6], 'EdgeColor', 'none'); 
  xlabel('Time (days)');
  ylabel('no. of infected');
  title('SIR Model - Infected Group');
  grid on;

% Recovered - Light Green
  subplot(3,1,3);
  area(t, R, 'FaceColor', [0.6 1 0.6], 'EdgeColor', 'none'); 
  xlabel('Time (days)');
  ylabel('no. of recovered');
  title('SIR Model - Recovered Group');
  grid on;

% ------------------------ PLOT : Complete SIR ---------------------------%
% Normalize the SIR data
  S_norm = S / N;
  I_norm = I / N;
  R_norm = R / N;

  figure;
  plot(t, S_norm, '-', 'Color', [0.1, 0.44, 0.9], 'LineWidth', 3);
  hold on;
  plot(t, I_norm, '--', 'Color', [0.95, 0.07, 0.09], 'LineWidth', 3);
  plot(t, R_norm, ':', 'Color', [0.37, 0.9, 0.18], 'LineWidth', 3);
  hold off;
% Create dummy variable to produce R_o value at the legend
  h_dummy = patch(nan, nan, 'k');
  
  xlabel('Time (days)');
  ylabel('Fraction of Population (density)');
  title('SIR Model Simulation (Normalized)'); 
  legend({'Susceptible', 'Infected', 'Recovered',...
          ['R_o = ' num2str(R_o, '%.2f')]}, 'Location', 'best');
  grid on;
  axis tight;
  set(gca, 'FontSize', 12); % Adjust font size
  box on;

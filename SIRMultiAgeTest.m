%-------------------------------------------------------------------------%
%          INITIALIZATION and PLOTS for SIR Age-Structured Model
%-------------------------------------------------------------------------%

% Initialize variables 
  N = 5000;                       % total number of particles
  NR = [0; 0; 0];                 % recovered particles per age group
  T_span = linspace(0, 89, 90);   % simulation duration time
  I_period = 7;                  % duration of infectious period
  infection_rate = 0.01;          % independent variable
  I0 = round(N * infection_rate); % initial infection rate
  X = [0.20; 0.63; 0.17];         % fraction of particles per age group
  
% Daily effective reproduction no. per age group for First Outbreak
  Re = [1.3; 0.85; 0.6];          % Kids, Adults, Elderly

% Run the SIR simulation
  [t, S, I, R] = SIRMultiAge(N, I0, NR, Re, X, T_span, I_period);

% -------------------------- RESULT ANALYSIS -----------------------------%

% [R G B] Color Palette 
  face_blue  = [0.6 0.8 1];
  face_red   = [1 0.6 0.6];
  face_green = [0.6 1 0.6];

% ------------------------ PLOT : Complete SIR ---------------------------%
  figure;
  plot(t, sum(S,2), '-', 'Color', [0.1, 0.44, 0.9], 'LineWidth', 3); 
  hold on;
  plot(t, sum(I,2), '-o', 'Color', [0.95, 0.07, 0.09], 'LineWidth', 2);
  plot(t, sum(R,2), '-.', 'Color', [0.37, 0.9, 0.18], 'LineWidth', 3);
  legend({'Susceptible','Infected','Recovered'}, 'FontSize', 12);
  xlabel('Time (days)');
  ylabel('Total Number of Individuals per Compartment');
  title('Disease Dynamics of COVID-19 Transmission (SIR)');
  grid on;
  hold off;
  box on;

% ------------------------ ANALYSIS: Peak Days ---------------------------%
% Find peak infections for each group
  [~, peak_kids] = max(I(:,1));
  [~, peak_adults] = max(I(:,2));
  [~, peak_elderly] = max(I(:,3));

% Convert peak day to real date
 % t = 1:100;
  start_date = datetime(2020,12,12);
  date_vector = start_date + days(t - 1);

  fprintf('\n📈 Peak Infection Days:\n');
  fprintf('  Kids:    Day %d (%s)\n', peak_kids...
                        , datestr(date_vector(peak_kids), 'mmm dd, yyyy'));
  fprintf('  Adults:  Day %d (%s)\n', peak_adults,...
                        datestr(date_vector(peak_adults), 'mmm dd, yyyy'));
  fprintf('  Elderly: Day %d (%s)\n', peak_elderly,...
                       datestr(date_vector(peak_elderly), 'mmm dd, yyyy'));

% Add peak markers on main plot
  figure(1); hold on;

  plot(peak_kids, sum(I(peak_kids,:)), 'o', ...
       'MarkerSize', 8,'MarkerEdgeColor', 'k',...
       'MarkerFaceColor', 'yellow','HandleVisibility','off');

  plot(peak_adults, sum(I(peak_adults,:)), 's', ...
      'MarkerSize', 8,'MarkerEdgeColor', 'k',...
      'MarkerFaceColor', 'yellow','HandleVisibility','off');

  plot(peak_elderly, sum(I(peak_elderly,:)), 'd', ...
      'MarkerSize', 8,'MarkerEdgeColor', 'k',...
      'MarkerFaceColor', 'yellow','HandleVisibility','off');
  legend('show');  

% --------------- PLOT: Individual Results per age group -----------------%
% Kids (0-11) years SIR dynamics
  figure;
  area(t, S(:,1), 'FaceColor', face_blue, 'EdgeColor', 'b',...
                  'FaceAlpha', 0.5); hold on;
  area(t, I(:,1), 'FaceColor', face_red, 'EdgeColor', 'r',...
                  'FaceAlpha', 0.5);
  area(t, R(:,1), 'FaceColor', face_green, 'EdgeColor', 'g',...
                  'FaceAlpha', 0.5);  

  xlabel('Time (days)');
  ylabel('Number of Individuals per Compartment');
  legend('S1','I1','R1');
  title('SIR Model for Kids (0-11) Age Groups');
  grid on;

% Adults (12-59) years SIR dynamics 
  figure;
  area(t, S(:,2), 'FaceColor', face_blue, 'EdgeColor', 'b',...
                  'FaceAlpha', 0.5); hold on;
  area(t, I(:,2), 'FaceColor', face_red, 'EdgeColor', 'r',...
                  'FaceAlpha', 0.5);
  area(t, R(:,2), 'FaceColor', face_green, 'EdgeColor', 'g',...
                  'FaceAlpha', 0.5);  

  xlabel('Time (days)');
  ylabel('Number of Individuals per Compartment');
  legend('S2','I2','R2');
  title('SIR Model for Adults (12-59) Age Groups');
  grid on;

% Elderly (60+) years SIR dynamics  
  figure;
  hold on;
  area(t, S(:,3), 'FaceColor', face_blue, 'EdgeColor', 'b',...
                  'FaceAlpha', 0.5); hold on;
  area(t, I(:,3), 'FaceColor', face_red, 'EdgeColor', 'r',...
                  'FaceAlpha', 0.5);
  area(t, R(:,3), 'FaceColor', face_green, 'EdgeColor', 'g',...
                  'FaceAlpha', 0.5);  

  xlabel('Time (days)');
  ylabel('Number of Individuals per Compartment');
  legend('S3','I3','R3'); 
  title('SIR Model for Adults (60+) Age Groups');
  grid on;

% --------------------- PLOT : Area Under the Curve ----------------------%
% Calculate the Area Under Curve for (infected) using Trapezoid rule 
  AUC_kids = trapz(t, I(:,1));    % Kids (0-11 years)
  AUC_adults = trapz(t, I(:,2));  % Adults (11-59 years)
  AUC_elderly = trapz(t, I(:,3)); % Elderly (60+ years)

% Display the labels for each age group
  fprintf('AUC for Kids (0-11): %.2f\n', AUC_kids);
  fprintf('AUC for Adults (12-59): %.2f\n', AUC_adults);
  fprintf('AUC for Elderly (60+): %.2f\n', AUC_elderly);

% Compare AUC of each age group and display the result
  if AUC_kids > AUC_adults && AUC_kids > AUC_elderly
     disp(['Kids are more likely to spread the infection based on' ...
           'AUC comparison.']);
  elseif AUC_adults > AUC_kids && AUC_adults > AUC_elderly
         disp(['Adults are more likely to spread the infection based'... 
               'on AUC comparison.']);
      else
         disp(['Elderly are more likely to spread the infection based' ...
               ' on AUC comparison.']);
  end

% Plot the AUC figure for infected population per age group
  figure;
  bar([AUC_kids, AUC_adults, AUC_elderly]);
  title('AUC (Cumulative Infected Population) for Different Age Groups');
  ylabel('Area Under the Curve (Infected Population)');
  xlabel('3 Distinct Age Groups');
  xticklabels({'Kids (0-11)', 'Adults (11-59)', 'Elderly (60+)'});
  grid on;


% ----------------- PLOT : Subplots by compartment -----------------------%
  figure;

% Susceptible
  subplot(3,1,1);
  plot(t, S(:,1), '-', 'Color', [0.1, 0.44, 0.9], 'LineWidth', 2);
  hold on;
  plot(t, S(:,2), '--', 'Color', [0.1, 0.44, 0.9], 'LineWidth', 2);
  plot(t, S(:,3), ':', 'Color', [0.1, 0.44, 0.9], 'LineWidth', 2);
  title('Susceptible by Age Group'); ylabel('No. of individuals');
  legend({'Kids','Adults','Elderly'}, 'FontSize', 10);
  grid on;

% Infected
  subplot(3,1,2);
  plot(t, I(:,1), '-', 'Color', [0.95, 0.07, 0.09], 'LineWidth', 2);
  hold on;
  plot(t, I(:,2), '--', 'Color', [0.95, 0.07, 0.09], 'LineWidth', 2);
  plot(t, I(:,3), ':', 'Color', [0.95, 0.07, 0.09], 'LineWidth', 2);
  title('Infected by Age Group'); ylabel('No. of individuals');
  legend({'Kids','Adults','Elderly'}, 'FontSize', 10);
  grid on;

% Recoeverd
  subplot(3,1,3);
  plot(t, R(:,1), '-', 'Color', [0.37, 0.9, 0.18], 'LineWidth', 2);
  hold on;
  plot(t, R(:,2), '--', 'Color', [0.37, 0.9, 0.18], 'LineWidth', 2);
  plot(t, R(:,3), ':', 'Color', [0.37, 0.9, 0.18], 'LineWidth', 2);
  title('Recovered by Age Group');...
        xlabel('Time (days)'); ylabel('No. of individuals');
  legend({'Kids','Adults','Elderly'}, 'FontSize', 10);
  grid on;

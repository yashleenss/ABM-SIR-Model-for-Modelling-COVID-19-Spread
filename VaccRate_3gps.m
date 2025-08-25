%-------------------------------------------------------------------------% 
%              Vaccination Percentages Three Age Groups 
%-------------------------------------------------------------------------%

% Load the Excel dataset for vaccination rate
  data = readtable('VaccinationRate.xlsx');

% Convert the 'VaccinationDate' to datetime format
  data.VaccinationDate = datetime(data.VaccinationDate,...
                             'InputFormat', 'dd/MM/yyyy');

% Time span 
  startDate = datetime('20-Dec-2020', 'InputFormat', 'dd-MMM-yyyy');
  endDate = datetime('26-Jul-2022', 'InputFormat', 'dd-MMM-yyyy');

% Three distinct age group categories based on vaccination rates
  kids = {'0-19'};
  adults = { '20-29', '30-39', '40-49', '50-59'};
  elderly = {'60-69', '70-79', '80+'};

% Percentage columns for plotting the graph
  doses = {'first_dose_percent', 'second_dose_percent',...
              'third_dose_percent', 'fourth_dose_percent'};
  titles = {'First Dose', 'Second Dose', 'Third Dose', 'Fourth Dose'};

% Ensure that all dose percentages are in the range of 0-100
  for j = 1:length(doses)

      % If values are above 1, assume they are in percent
        if max(data.(doses{j})) > 1  
           data.(doses{j}) = data.(doses{j});  % keep unchanged
        else
           data.(doses{j}) = data.(doses{j}) * 100;  % convert to percent
        end
   end

% To aggregate age group data
  aggregateAgeGroup = @(groupNames) varfun(@mean,...
                       data(ismember(data.age_group, groupNames), :),...
                       'InputVariables', doses, 'GroupingVariables',...
                       'VaccinationDate');

% Aggregate data for the 3 age groups
  kidsData = aggregateAgeGroup(kids);
  adultsData = aggregateAgeGroup(adults);
  elderlyData = aggregateAgeGroup(elderly);

% Define key dates for 1st, 2nd, and 3rd vaccination campaign
  majorDates = datetime({'08-Jan-2021', '13-Mar-2021', '01-Aug-2021'}); 

% Initialize figure
  figure;

  for j = 1:length(doses)
      subplot(2, 2, j);
      hold on;

     % Plot combined age groups
       plot(kidsData.VaccinationDate, kidsData.(sprintf('mean_%s',...
                      doses{j})), 'b-', 'DisplayName', 'Kids (0-19)');
       plot(adultsData.VaccinationDate, adultsData.(sprintf('mean_%s',...
                   doses{j})), 'r-', 'DisplayName', 'Adults (11-59)');
       plot(elderlyData.VaccinationDate, elderlyData.(sprintf('mean_%s',...
                      doses{j})), 'g-', 'DisplayName', 'Elderly (60+)');

      % Add vertical lines at major dates
        for k = 1:length(majorDates)
            xline(majorDates(k), '--k', 'HandleVisibility', 'off');
        end

      % X-axis ticks and labels
        xticks([datetime('20-Dec-2020'), majorDates]); 
        xticklabels({'Dec 20, 2020','Jan 8, 2021', 'Mar 13', 'Aug 1'});
        xlabel('Date');
        ylabel('% of Population Vaccinated');
        title(titles{j});
        ylim([0 100]);
        xlim([startDate, endDate]);
        datetick('x', 'dd-mmm-yyyy', 'keeplimits');
        legend('show', 'Location', 'northeast');
        grid on;
        hold off;
  end

% Title for the entire figure
  sgtitle('Vaccination Dose Percentages Over Time by Age Groups');

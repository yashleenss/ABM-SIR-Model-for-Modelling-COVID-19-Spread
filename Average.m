function Avg = Average(P_H0)

  N = 5000;
% Age Group distribution (percentages)
  elderly_percentage = 0.17;     % 17% are elderly
  adult_percentage = 0.63;       % 63% are adults
  kids_percentage = 0.20;        % 20% are kids

% Number of individuals in each group
  n_elderly = round(N * elderly_percentage);
  n_adults = round(N * adult_percentage);
  n_kids = N - n_elderly - n_adults; % Remaining are kids
    
% Separate the P_H vector into 3 age groups
  kids = P_H0(1:n_kids, :); 
  adults = P_H0(n_kids+1:n_kids+n_adults, :); 
  elderly = P_H0(n_kids+n_adults+1:end, :); 
    
% Calculate the average summation per day for each group
  Avg.kids = sum(kids, 1) / n_kids;
  Avg.adults = sum(adults, 1) / n_adults;
  Avg.elderly = sum(elderly, 1) / n_elderly;

end
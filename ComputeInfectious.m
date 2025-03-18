function infection = ComputeInfectious(tI, t, variant)

% Variant type and calculate the infectious period
  switch variant
  case 'Alpha and Delta'

    % Use element-wise logical operators RHS and LHS are vectors
    % Infection period for an individual is between 4th and 7th day
      if ((t-tI) >= 4) && ((t-tI) <= 7)
          infection = 1;
      else
          infection = 0;
      end 
        
 
  % Infection period for an individual is between 2nd and 5th day
    case 'Omicron'
    if ((t-tI) >= 2) && ((t-tI) <= 5)
          infection = 1;
      else
          infection = 0;
    end  

    otherwise
          infection = 0;
  end
end

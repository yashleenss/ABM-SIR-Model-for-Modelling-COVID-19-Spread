%-----------------------------------------------------%
% Calculating the initial infection period and variant
%-----------------------------------------------------%

function infection = ComputeInfectious(tI, t, variant)

      % Infectious period according to variant type
        switch variant

      % Use element-wise logical operators RHS and LHS are vectors
      % Infection period for an individual is between 4th and 7th day
        case 'Alpha and Delta'
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

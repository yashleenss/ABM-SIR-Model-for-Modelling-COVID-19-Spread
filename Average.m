%-------------------------------------------------------------------------%
%       Calculating AVERAGE infection probability per age group 
%-------------------------------------------------------------------------%

function Avg = Average(age_group, P_H0)
    
        % Calculate the average summation per day for each group
          Avg.elderly = mean(P_H0(age_group == 1, :), 1);
          Avg.adults = mean(P_H0(age_group == 2, :), 1);
          Avg.kids = mean(P_H0(age_group == 3, :), 1);

end

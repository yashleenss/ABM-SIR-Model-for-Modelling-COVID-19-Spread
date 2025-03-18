function VE = VaccineEffectiveness(Tv, t)

% t1: denotes the first day since infection started (106)
% Tv: indicates the day of immunization/vaccination for each particle
% Initialize VEdays and VElevel

  VElevel = 0.9; % Maximum VElevel is 90%
  VEdays = [7, 28, 150, 180]; % Days after which VE changes
   
% First, check if the particle is unvaccinated (Tv should be infinity)
 if isinf(Tv)
     VE = 0; % If unvaccinated, set VE=0
     return;
 end

% Compute vaccine effectiveness based on the days since vaccination
% (t-Tv) will compute the number of days since vaccination
    
  if (t-Tv) <= VEdays(1)
     VE = 0; % No effectiveness before first 7 days
  else
      if isscalar(t-Tv) % If it's a single value, use scalar logical operators
         if (t-Tv) > VEdays(1) && (t-Tv) <= VEdays(2) 
             % Linear increase from 0% to 90% over the next 21 days   
               VE = VElevel * ((t-Tv - VEdays(1)) / (VEdays(2) - VEdays(1)));
         else
             if (t-Tv) > VEdays(2) && (t-Tv) <= VEdays(3)
                % VE remains at 90% until 150 days after 1 week of second dose
                  VE = VElevel; % Maximum effectiveness between 28 and 150 days
             else
                 if (t-Tv) > VEdays(3) && (t-Tv) <= VEdays(4)
                    % Linear decrease from 90% to 30% over the next 180 days
                      VE = VElevel * (1 - ((t-Tv - VEdays(3)) / (3 * (VEdays(4) - VEdays(3)))));
                 else
                     % Decline in VE after 180 days (should be 30% or/6)
                       VE = 2 * (VElevel / 3); 
                  end
              end
          end
      end 
  end
end
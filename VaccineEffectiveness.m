%-------------------------------------------------------------------------%
%    Calculating Vaccine Effectiveness based on vaccination/recovery
%-------------------------------------------------------------------------%

function VE = VaccineEffectiveness(Ti, t, VElevel, VEdays, Immunity)

% t1 = denotes the first day since infection started since Aug/28/2020
% Ti = indicates the day of immunization/vaccination for each particle
% (t-Ti) will compute the number of days since vaccination/immunization
   
% First, check if the particle is unvaccinated (Tv should be infinity)
  if Ti == Inf
     VE = 0;   % if unvaccinated, set VE=0
     return;
  end

%------------------------------ IMMUNITY MODEL ---------------------------%

  if Immunity == 1

    % Compute VE decay based on natural infection (decay model)
     if (t-Ti) < VEdays(2)  
        VE = VElevel;

     elseif (t-Ti) > VEdays(2) && (t-Ti) <= VEdays(4)
         
            % Linear decrease from 0.9% to 0.3% with -0.6/152 slope
              m = (0.3 - VElevel) / (VEdays(4) - VEdays(2));
              VE = m * ((t - Ti) - VEdays(2)) + VElevel;    
     else
            % Maintain ~60% after 180 days for Alpha and Delta variant
            % maintain ~30% after 180 days for Omicron variant
              VE = (VElevel / 3);  

            
     end

%---------------------------- VACCINATION MODEL --------------------------% 
  else

% Compute VE growth based on vaccination (growth model)
  if (t-Ti) <= VEdays(1)
     VE = 0; % No effectiveness before first 7 days

  else
      if isscalar(t-Ti) % for a single value, use scalar logical operators
         if (t-Ti) > VEdays(1) && (t-Ti) <= VEdays(2) 

         % Linear increase from 0% to 90% over the next 21 days   
           VE = VElevel * (((t-Ti) - VEdays(1)) / (VEdays(3) - VEdays(2)));

         else
             if (t-Ti) > VEdays(2) && (t-Ti) <= VEdays(3)

               % VE remains at 90% for 150 days after 1 week of second dose
               % Maximum effectiveness between 28 and 150 days
                 VE = VElevel;

             else
                 if (t-Ti) > VEdays(3) && (t-Ti) <= VEdays(4)

                   % Linear decrease from 90% to 30% over the next 180 days
                     VE = VElevel * (1 - (((t-Ti) - VEdays(3)) /...
                           (3 * (VEdays(4) - VEdays(3)))));
                 else
               
                 % Maintain ~60% after 180 days for Alpha and Delta variant
                 % maintain ~30% after 180 days for Omicron variant
                     VE = (VElevel / 3); 

                  end
              end
          end
      end 
  end
  end
end 


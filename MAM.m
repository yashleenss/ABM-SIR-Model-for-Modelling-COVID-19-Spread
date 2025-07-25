%-------------------------------------------------------------------------%
%         MAM Algorithm for calculating the infection probability
%-------------------------------------------------------------------------%

  function results = MAM(N, L, H, X, Y, T_span, Tv, Ti, tI, VE, VElevel,...
                     VEdays, VE_fixed, ngps, age_group,...
                     infection, radius, RecoveryTime, Immunity,...
                     ImmunityTime, daily_new_elderly, R_t_vector,...
                     daily_new_adults, daily_new_kids, P_H_tspan)

           % Initialize arrays to store the no. of susceptible and infected
             daily_susceptible_idx = zeros(1,T_span);
             daily_infected_idx = zeros(1,T_span);
             daily_recovery_idx = zeros (1, T_span);

           % Initialize storage vectors for time-evolving probabilities 
             P_H0 = zeros(N, T_span);
             VE_all = zeros(1, T_span);       % Avg VE for N population
             VE_infected = zeros(1, T_span);  % Avg VE among infected only

  % Loop through each day to update positions and check infections 
    for t = 1:T_span

        % Updating the variant type
          variant = variant_type(t);

        % Periodic Boundary length(Lo) of the simulation area
          Lo = L * sqrt(R_t_vector(t) / 1.25);
          sigma = Lo / 2;       % Standard Deviation
          sigma_r = sigma/1;    % for the infection probability calculation


        %---------------------- PARTICLE POSITION ------------------------%
        % Compute daily displacements with Brownian motion
          deltaX = normrnd(0, sigma, [N, 1]);
          deltaY = normrnd(0, sigma, [N, 1]);

        % Wrapped around boundary using modulo (no particle is leaving Lo)
          X = mod(X + deltaX, L);
          Y = mod(Y + deltaY, L);

        % Store temporary probabilities for the current time step
          P_H = zeros(N, 1);     % initialization for probability

        %------------------------- RECOVERY STATE ------------------------%
          for j = 1:N
              if H(j) == 1 
                 if tI(j)+ RecoveryTime(j) <= t
                   H(j) = 2;         % return to recovered state
                   Immunity(j) = 1;  % recovery based immunity
                   tI(j) = 0;        % no infection day
                   infection(j) = 0; % no longer actively infectious
                   Ti(j) = t;        % reset recovery day to calculate VE
                 end
              end
          end

        %-------------------Update Recovered Individuals------------------%
        % Reset these particles to susceptible after the temporary immunity

          for j = 1:N
              if H(j) == 2 && (t - Ti(j)) >= ImmunityTime(j)  
                  H(j) = 0;      
                  infection(j) = 0;
              end
          end

        %--------------------- VACCINATION DOSE CHECK --------------------%
          for dose = 1:3

             % update any new dose recieved at current time t
               newly_vaccinated = find (Tv(:, dose) == t); 
               for j = newly_vaccinated
                   if Immunity(j) ~= 1  % check for no recovery-immunity
                      Ti(j) = t;        % reset recovery day 
                      Immunity(j) = 0;  % mark as vaccine-immunity
                   end
               end
          end

%------------------------ INFECTION PROBABILITY CHECK --------------------%
% Begin by checking the health status for each particle

   for i = 1:N

     % Check if particle i is vaccinated and update VE
       VE(i) = VaccineEffectiveness(Ti(i), t, VElevel,VEdays, Immunity(i));

     % If the particle is healthy (H=0)
       if H(i) == 0 
           P_H(i) = 0;     % default infection probability
   
        % Check the infection probability for the next particle
           for j = 1:N

             % If the other particle is infected (H = 1)
               if i~=j 
                  if H(j) == 1 && tI(j) < t
                 
              % Calculate the distance between i and j particles
                dx = abs(X(i) - X(j));
                dy = abs(Y(i) - Y(j));
                  
              % Shortest wrapped around distance across the PB  
                dx = min(dx, L - dx);
                dy = min(dy, L - dy);
                distance = dx^2 + dy^2; 
                
                   % Check the infection radius 
                     if distance <= radius
                        exp_term = exp(-distance / (2 * (sigma_r)^2));
                    
                     % Infection day for each infected particle
                       infection(j) = ComputeInfectious(tI(j),t, variant);

                      %---------------------------------------------------% 

                      % Calculating the probability of getting an infection
                       P_H(i)=P_H(i)+(infection(j)*(1 - VE(i))*(exp_term));

                      %---------------------------------------------------%  
             
                     end
                  end 
               end 
           end
                    
           % Ensure Clamping values within [0, 1]
             P_H(i) = max(0, min(1, P_H(i))); 

           % Uniformly Distributed number within [0,1] for each particle 
             xi = rand;
                
           % Only update for the NEWLY infected individuals
             if xi < P_H(i) && H(i) == 0 
                 H(i) = 1;         % mark these particle as infected            
                 tI(i) = t;        % reset infection day 
                 VE(i) = 0;        % reset no vaccine effectiveness
                 infection(i) = 1; % reset the particle to be infectious
             end
       end
   end
           
  % Store the daily count /age group (susceptible, infected & recovered)
    daily_susceptible_logical = (H == 0); 
    daily_susceptible_idx(t) = sum (daily_susceptible_logical);

    daily_infected_logical = (H == 1);   % (N*1) array vector
    daily_infected_idx(t) = sum(daily_infected_logical);  % scalar

    daily_new_elderly(t) = sum(daily_infected_logical & (age_group == 1));
    daily_new_adults(t) = sum(daily_infected_logical & (age_group == 2));
    daily_new_kids(t) = sum(daily_infected_logical & (age_group == 3));

    daily_recovery_logical = (H == 2);
    daily_recovery_idx(t) = sum(daily_recovery_logical);

  % Calling the avergae function
    P_H0(:,t) = P_H;
    Avg = Average(age_group, P_H0);

 % Separate VE groups (1,2,3) and store their H values
   VE_group = discretize(VE, [-Inf, VE_fixed, Inf]);
   for i = 1:ngps

       % fraction of individuals in VE group i who are infected 1
         P_H_tspan(i,t) = mean(H(VE_group==i)==1);
   end

 % Avergae vaccine effectiveness over all N individuals at time t
   VE_all(t) = mean(VE);
   
 % Average vaccine effectiveness only for infected individuals at time t
    if any(daily_infected_logical)
       VE_infected(t) = mean(VE(daily_infected_logical));
    else
       VE_infected(t)=0;
    end

     end % end time span loop 

  % Store the results from each unique run
    results.H = H;
    results.P_H = P_H;
    results.Avg = Avg;
    results.VE_all = VE_all;
    results.P_H_tspan = P_H_tspan;
    results.VE_infected = VE_infected;
    results.daily_new_kids = daily_new_kids;
    results.daily_new_adults = daily_new_adults;
    results.daily_new_elderly = daily_new_elderly;
    results.daily_infected_idx = daily_infected_idx;
    results.daily_recovery_idx = daily_recovery_idx;
    results.daily_susceptible_idx = daily_susceptible_idx;
  end 

  % Select a variant type to based on the time span
    function variant = variant_type(t)
              if t >= 450
                    variant = 'Omicron';
              else
                    variant = 'Alpha and Delta';
              end
      end

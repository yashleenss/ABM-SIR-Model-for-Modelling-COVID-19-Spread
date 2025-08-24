%-------------------------------------------------------------------------%
%         Calculating the initial Health Status of each particle
%-------------------------------------------------------------------------%

 function [H] = HealthStatus(N, NRc)

    % N = Total number of particles
    % NSI = Number of initially susceptible particles
    % NI = Number of initially infected particles
    % NRc = Number of recovered particles
   
    % Initialization
      infection_rate = 0.01; % independent variable
      NI = round(N * infection_rate); 
      NSI = N - NRc;
      H = zeros(N, 1);      % Health status (0: susceptible, 1; infected)

    % Choose NI particles to be initially infected from the first NSI
      infected_particles = randperm(NSI, NI);  % randomly select particles
      
    % Set health status of selected particles to be infected
      H(infected_particles) = 1;
      
    % The remaining particles are susceptible, H is initialized with zeros.
    % No need to explicitly set them to 0.
    
 end

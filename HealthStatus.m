function [X, Y, H] = HealthStatus(N, NI, NRc, Lo)
    % N: Total number of particles
    % NI: Number of initially infected particles
    % NRc: Number of recovered particles
    % Lo: The length of one side of the boundary (particles randomly placed within 0 to Lo)
    
    % Total number of susceptible particles (excluding those recovered)
      Ns = N - NRc;

    % Initialize arrays to store the coordinates and health status
      X = Lo * rand(N, 1);  % X-coordinates, random between 0 and Lo
      Y = Lo * rand(N, 1);  % Y-coordinates, random between 0 and Lo
      H = zeros(N, 1);      % Health status (0 for susceptible, 1 for infected)

    % Choose N1 particles to be initially infected from the first Ns1
      infected_particles = randperm(Ns, NI);  % Randomly select N1 pasticles within Ns1

    % Set health status of selected particles to infected
      H(infected_particles) = 1;

    % The remaining particles are susceptible
    % No need to explicitly set them to 0 because H is initialized with zeros

end
%-------------------------------------------------------------------------%
%               SIR Compartmental Age-Structured Model
%-------------------------------------------------------------------------%

function [t, S, I, R] = SIRMultiAge(N, I0, NRC, Re, X, T_span, I_period)

         % Number of distinct groups (age groups)
           n = numel(X);
    
         % Define initial conditions for each age group
           NI = I0 .* X ;       % Initial infected per group
           NRC = NRC .* X ;       % Initial recovered per group
           NSI = N - NI - NRC;    % Initial susceptible per group
           Y0 = [NSI; NI; NRC];   % Initial conditions t solve dydt
     
         % Define the differential equations
           function dydt = SIR_eq(~, Y)
                       n = numel(Y) / 3;
                       S = Y(1:n);
                       I = Y(n+1:2*n);
                       R = Y(2*n+1:end);
        
         % Calculate infection and recovery rate for each age group
           mu =  (1 / I_period);
           beta = Re * (mu / N);

         % Define differential equations for each n-th age group
           dSdt = -beta .* S .* sum(I);
           dIdt = - dSdt - (mu * I);
           dRdt = mu * I;
        
           dydt = [dSdt; dIdt; dRdt];
           end

         % Solving ODE (ensuring non-negative solutions)
           options = odeset('NonNegative', 1:n*3); 
           [t, Y] = ode45(@SIR_eq, T_span, Y0, options);
    
         % Separate the results into S, I, R vectors 
           S = Y(:, 1:n);
           I = Y(:, n+1:2*n);
           R = Y(:, 2*n+1:end);
end

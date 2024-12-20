%% Part 1: modeling disease spread using SIR model
format shortG; clear; clc

% Initializing Variable Definitions
St1 = zeros(1,1);
It1 = zeros(1,1);
Rt1 = zeros(1,1);
St2 = zeros(1,1);
It2 = zeros(1,1);
Rt2 = zeros(1,1);
St3 = zeros(1,1);
It3 = zeros(1,1);
Rt3 = zeros(1,1);
St(1,1) = 990; % Number of susceptible individuals at time t
It(1,1) = 10; % Number of infected individuals at time t
Rt(1,1) = 0; % Number of recovered individuals at time t

N = St + It + Rt; % Total Population
h = 1; % Step size in days
T = 100; % Total simulation time in days (0 -> 100)
days = 1:h:100;

for t = 1:h:T % For loop calculation for the standard flu 

    % initially given variables
    Beta = 0.3; %Transmission Rate
    Gamma = 0.1; %Recovery Rate

    % Set of nonlinear differential equations defined
    dSdt = @(t, St, It, Rt) -(Beta / N) .* St .* It;
    dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
    dRdt = @(t, St, It, Rt) Gamma .* It;

    % Runge Kutta Calculations in which each diff eq is evaluated  
    k1S = dSdt(t, St(t,1), It(t,1), Rt(t,1));
    k1I = dIdt(t, St(t,1), It(t,1), Rt(t,1));
    k1R = dRdt(t, St(t,1), It(t,1), Rt(t,1));

    k2S = dSdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);
    k2I = dIdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);
    k2R = dRdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);

    k3S = dSdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);
    k3I = dIdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);
    k3R = dRdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);

    k4S = dSdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);
    k4I = dIdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);
    k4R = dRdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);


    St(t + 1,1) = St(t,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t + 1,1) = It(t,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t + 1,1) = Rt(t,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

    % Values for susceptible, infected, and recovered persons are tabulated and synthesized into a vector

    St1(t,1) = St(t,1);
    It1(t,1) = It(t,1);
    Rt1(t,1) = Rt(t,1);

end


for t = 1:h:T

    Beta = 1;
    Gamma = 0.1;
    
    dSdt = @(t, St, It, Rt) -(Beta / N) .* St .* It;
    dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
    dRdt = @(t, St, It, Rt) Gamma .* It;

    k1S = dSdt(t, St(t,1), It(t,1), Rt(t,1));
    k1I = dIdt(t, St(t,1), It(t,1), Rt(t,1));
    k1R = dRdt(t, St(t,1), It(t,1), Rt(t,1));

    k2S = dSdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);
    k2I = dIdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);
    k2R = dRdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);

    k3S = dSdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);
    k3I = dIdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);
    k3R = dRdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);

    k4S = dSdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);
    k4I = dIdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);
    k4R = dRdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);


    St(t + 1,1) = St(t,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t + 1,1) = It(t,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t + 1,1) = Rt(t,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

    St2(t,1) = St(t,1);
    It2(t,1) = It(t,1);
    Rt2(t,1) = Rt(t,1);

end

% Begin 
for t = 1:h:T

    % Re-initialize givens
    Beta = 2;
    Gamma = 0.2;
    
    % Next set of nonlinear differential equations defined
    dSdt = @(t, St, It, Rt) -(Beta / N) .* St .* It;
    dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
    dRdt = @(t, St, It, Rt) Gamma .* It;
    
    % Define all constants for the RK method, for each model
    k1S = dSdt(t, St(t,1), It(t,1), Rt(t,1));
    k1I = dIdt(t, St(t,1), It(t,1), Rt(t,1));
    k1R = dRdt(t, St(t,1), It(t,1), Rt(t,1));

    k2S = dSdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);
    k2I = dIdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);
    k2R = dRdt(t + 0.5 * h, St(t,1) + 0.5 * k1S * h, It(t,1) + 0.5 * k1I * h, Rt(t,1) + 0.5 * k1R * h);

    k3S = dSdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);
    k3I = dIdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);
    k3R = dRdt(t + 0.5 * h, St(t,1) + 0.5 * k2S * h, It(t,1) + 0.5 * k2I * h, Rt(t,1) + 0.5 * k2R * h);

    k4S = dSdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);
    k4I = dIdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);
    k4R = dRdt(t + h, St(t,1) + k3S * h, It(t,1) + k3I, Rt(t,1) + k3R);

    % RK method, employing constants created ealrier
    St(t + 1,1) = St(t,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t + 1,1) = It(t,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t + 1,1) = Rt(t,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

    St3(t,1) = St(t,1);
    It3(t,1) = It(t,1);
    Rt3(t,1) = Rt(t,1);

end


figure(1) % Plot Seasonal Influenza
hold on
grid on

plot(days,St1,'red','LineWidth',2)
plot(days,It1,'green','LineWidth',2)
plot(days,Rt1,'blue','linewidth',2)

title('Seasonal Influenza, ß = 0.3 y = 0.1')
xlabel('Time in Days')
ylabel('People')
legend('Susceptible','Infected','Recovered')
hold off


figure(2) % Plot COVID-19
hold on
grid on

plot(days,St2,'r','LineWidth',2)
plot(days,It2,'g','LineWidth',2)
plot(days,Rt2,'b','linewidth',2)

title('COVID-19, ß = 1 y = 0.1')
xlabel('Time in Days')
ylabel('People')
legend('Susceptible','Infected','Recovered')
hold off


figure(3) % Plot Measles
hold on
grid on

plot(days,St3,'r','LineWidth',2)
plot(days,It3,'g','LineWidth',2)
plot(days,Rt3,'b','linewidth',2)

title('Measles, ß = 2 y = .2')
xlabel('Time in Days')
ylabel('People')
legend('Susceptible','Infected','Recovered')
hold off


%% Discussion Section:

% The beta value is the transmission rate, or how often susceptible
% individuals get infected while gamma is the recovery rate, or how often the
% affected individuals recover. When the beta parameter increases, the
% number of susceptible individuals decreases significantly quicker, the
% overall recovery rate is considerably quicker, and the infected
% population spikes to a greater amplitude earlier. When the gamma parameter
% doubles along with the beta parameter, the susceptible, infected, and
% recovered populations all have similar behavior, except they all exhibit
% the end behavior twice as quickly. the derivatives are much sharper
% (change in the populations occurs much quicker). 

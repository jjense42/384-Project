%% Part 2: Linear/Quadratic Interpolation
format shortG; clear; clc

% Variable Definitions
St1 = zeros(1,1);
It1 = zeros(1,1);
Rt1 = zeros(1,1);
St(1,1) = 990; % Number of susceptible individuals at time t
It(1,1) = 10; % Number of infected individuals at time t
Rt(1,1) = 0; % Number of recovered individuals at time t

N = St + It + Rt; % Total Combined Population
h = 1; % Step size in days
T = 100; % Total simulation time in days (0 -> 100)
days = 1:h:100;

for t = 1:h:T
    % Define 1st set of givens
    Beta = 0.3; % Transmission Rate
    Gamma = 0.1; % Recovery Rate

    % Change in individual populations
    dSdt = @(t, St, It, Rt) - (Beta / N) .* St .* It; 
    dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
    dRdt = @(t, St, It, Rt) Gamma .* It;

    % Initiate constants for each population individually
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

    % RK Method, employing constants previously defined
    St(t + 1,1) = St(t,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t + 1,1) = It(t,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t + 1,1) = Rt(t,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

    St1(t,1) = St(t,1);
    It1(t,1) = It(t,1);
    Rt1(t,1) = Rt(t,1);

end

St1True(:,:) = St1(:,:);
It1True(:,:) = It1(:,:);
Rt1True(:,:) = Rt1(:,:);

for t = 1:h:T
    if mod(t,2) ~= 0
        St1(t,1) = 0;
        It1(t,1) = 0;
        Rt1(t,1) = 0;
    end
end

% Synthesize all linear values for each model into matrix
St1linear(:,:) = St1(:,:);
It1linear(:,:) = It1(:,:);
Rt1linear(:,:) = Rt1(:,:);
% Synthesize all quadratic values for each model into matrix
St1Quad(:,:) = St1(:,:);
It1Quad(:,:) = It1(:,:);
Rt1Quad(:,:) = Rt1(:,:);

for t = 3:h:T - 1
    if mod(t,2) ~= 0
        St1linear(t,1) = ((t - (t + 1))/((t - 1) - (t + 1))) * St1(t - 1,1) ...
        + ((t - (t - 1))/((t + 1) - (t - 1))) * St1(t + 1,1);

        It1linear(t,1) = ((t - (t + 1))/((t - 1) - (t + 1))) * It1(t - 1,1) ...
        + ((t - (t - 1))/((t + 1) - (t - 1))) * It1(t + 1,1);
        
        Rt1linear(t,1) = ((t - (t + 1))/((t - 1) - (t + 1))) * Rt1(t - 1,1) ...
        + ((t - (t - 1))/((t + 1) - (t - 1))) * Rt1(t + 1,1);
    end
    
    St1linear(1,1) = (-3)/(-2) * St1(2,1) + (-1)/(2) * St1(4,1);
    It1linear(1,1) = (-3)/(-2) * It1(2,1) + (-1)/(2) * It1(4,1);
    Rt1linear(1,1) = (-3)/(-2) * Rt1(2,1) + (-1)/(2) * Rt1(4,1);

end

% Begin Lagrange quadratic interpolation
for t = 3:h:T - 2
    if mod(t,2) ~= 0
        St1Quad(t,1) = ((t - (t + 1)) * (t - (t + 3)) / (((t - 1) - (t + 1)) * ((t - 1) - (t + 3)))) * St1(t - 1,1)...
        + ((t - (t - 1)) * (t - (t + 3)) / (((t + 1) - (t - 1)) * ((t + 1) - (t + 3)))) * St1(t + 1,1)...
        + ((t - (t - 1)) * (t - (t + 1)) / (((t + 3) - (t - 1)) * ((t + 3) - (t + 1)))) * St1(t + 3,1);

        It1Quad(t,1) = ((t - (t + 1)) * (t - (t + 3)) / (((t - 1) - (t + 1)) * ((t - 1) - (t + 3)))) * It1(t - 1,1)...
        + ((t - (t - 1)) * (t - (t + 3)) / (((t + 1) - (t - 1)) * ((t + 1) - (t + 3)))) * It1(t + 1,1)...
        + ((t - (t - 1)) * (t - (t + 1)) / (((t + 3) - (t - 1)) * ((t + 3) - (t + 1)))) * It1(t + 3,1);

        Rt1Quad(t,1) = ((t - (t + 1)) * (t - (t + 3)) / (((t - 1) - (t + 1)) * ((t - 1) - (t + 3)))) * Rt1(t - 1,1)...
        + ((t - (t - 1)) * (t - (t + 3)) / (((t + 1) - (t - 1)) * ((t + 1) - (t + 3)))) * Rt1(t + 1,1)...
        + ((t - (t - 1)) * (t - (t + 1)) / (((t + 3) - (t - 1)) * ((t + 3) - (t + 1)))) * Rt1(t + 3,1);

    end

    % Interpolating boundary (right/left ending) values of each model manually
    St1Quad(1,1) = (15/8) * St1(2,1) + (-5/4) * St1(4,1) + (3/8) * St1(6,1);
    St1Quad(99,1) = (-1/8) * St1(96,1) + (3/4) * St1(98,1) + (3/8) * St1(100,1);
    It1Quad(1,1) = (15/8) * It1(2,1) + (-5/4) * It1(4,1) + (3/8) * It1(6,1);
    It1Quad(99,1) = (-1/8) * It1(96,1) + (3/4) * It1(98,1) + (3/8) * It1(100,1);
    Rt1Quad(1,1) = (15/8) * Rt1(2,1) + (-5/4) * Rt1(4,1) + (3/8) * Rt1(6,1);
    Rt1Quad(99,1) = (-1/8) * Rt1(96,1) + (3/4) * Rt1(98,1) + (3/8) * Rt1(100,1);

end

% Error calculations
El2_Sl = sqrt(sum(St1linear - St1True).^2 ./N);
El2_Il = sqrt(sum(It1linear - It1True).^2 ./N);
El2_Rl = sqrt(sum(Rt1linear - Rt1True).^2 ./N);
% Combine linear error values into a matrix
linerr = [El2_Sl; El2_Il; El2_Rl];

El2_Sq = sqrt(sum(St1Quad - St1True).^2 ./N);
El2_Iq = sqrt(sum(It1Quad - It1True).^2 ./N);
El2_Rq = sqrt(sum(Rt1Quad - Rt1True).^2 ./N);
% Combine quadratic error values into a matrix
quaderr = [El2_Sq; El2_Iq;El2_Rq];

% Creating a table to compile all error results
ErrorTable = table([linerr],[quaderr],'VariableNames', {'Linear Error', 'Quadratic Error'});
disp(ErrorTable);

%% Discussion section: 
% The quadratic Lagrangian method produces a significantly smaller 
% error than the linear Lagrangian method. The quadratic method provides a 
% better fit for functions with notable curvature, such as the specific 
% functions we have approximated in this part. Both of the methods used the
% same amount of data points, so the effect of data amounts and step sizes 
% were neglected to compare the accuracy of both methods specifically. 

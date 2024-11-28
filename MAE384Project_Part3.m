%% Part 3: Least Squares
format shortG; clear; clc

% Initializing Given Variables
St1 = zeros(1,1);
It1 = zeros(1,1);
Rt1 = zeros(1,1);

St(1,1) = 990; % Number of susceptible individuals at time t
It(1,1) = 10; % Number of infected individuals at time t
Rt(1,1) = 0; % Number of recovered individuals at time t

N = St + It + Rt; % Total Population
h = 1; % Step size in days
T = 30; % Total simulation time in days (0 -> 100)
days = 1:h:100;


for t = 1:h:T
    
    % Defining additional given values
    Beta = 0.3; %Transmission Rate
    Gamma = 0.1; %Recovery Rate

    % Defining diff eqs for Susceptible, Infected, and Recovered populations
    dSdt = @(t, St, It, Rt) -(Beta / N) .* St .* It;
    dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
    dRdt = @(t, St, It, Rt) Gamma .* It;

    % Defining all needed constants
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

    % Approximating, using all constants produced previously
    St(t + 1,1) = St(t,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t + 1,1) = It(t,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t + 1,1) = Rt(t,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

    St1(t,1) = St(t,1);
    It1(t,1) = It(t,1);
    Rt1(t,1) = Rt(t,1);

end

% linear least squares with T = 30
T = 30;
x=days(1:T);
y=log((It1(1:T))');

% Solve for ln I(t)=ln I(0)+kt. A1=k A0=ln I(0)
A1= (T*sum(x.*y)-(sum(x)*sum(y)))/(T*sum(x.^2)-(sum(x)^2));
A0= (sum(y)/T)-((A1/T)*(sum(x)));
I0=exp(A0);
b=(A1+Gamma)*(N/St(1,1));

% Linear least squares with T = 10
T = 10;
x1=days(1:T);
y1=log((It1(1:T))');

% Solving for ln I(t)=ln I(0)+kt. A11=k A01=ln I(0)
A11= (T*sum(x1.*y1)-(sum(x1)*sum(y1)))/(T*sum(x1.^2)-((sum(x1))^2));
A01= ((sum(y1))/T)-(A11*(sum(x1)/T));
I01=exp(A01);
b1=(A11+Gamma)*(N/St(1,1));

% Display all results
disp('Estimimated ß for T = 30:');
disp(b);
disp('Estimimated I(O) for T = 30:');
disp(I0);
disp('Estimimated ß for T = 10:');
disp(b1);
disp('Estimimated I(O) for T = 10:');
disp(I01);

%% Discussion Section

% As you decrease the T value (number of days), the estimated values for 
% I(0) and beta converge to the true values. Reducing the number of data 
% points, the I(t) curve approaches an exponential model within those boundaries.
% Therefore, the exponential approximation will  predict the actual data points 
% more accurately. 

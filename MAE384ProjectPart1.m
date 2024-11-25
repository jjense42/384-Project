%Runge-Kutta Method
% i can edit code here
St(1,1) = 990; %Number of susceptible individuals at time t
It(1,1) = 10; %Number of infected individuals at time t
Rt(1,1) = 0; %Number of recovered individuals at time t
N = St + It + Rt; %Total Population
Beta = 0.3; %Transmission Rate
Gamma = 0.1; %Recovery Rate
h = 1; %Step size in days
T = 100; %Total simulation time in days (0 -> 100)

dSdt = @(t, St, It, Rt) -(Beta / N) .* St .* It;
dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
dRdt = @(t, St, It, Rt) Gamma .* It;

for t = 1:h:T
    
    k1S = dSdt(t, St, It, Rt);
    k1I = dIdt(t, St, It, Rt);
    k1R = dRdt(t, St, It, Rt);

    k2S = dSdt(t + 0.5 * h, St + 0.5 * k1S * h, It + 0.5 * k1I * h, Rt + 0.5 * k1R * h);
    k2I = dIdt(t + 0.5 * h, St + 0.5 * k1S * h, It + 0.5 * k1I * h, Rt + 0.5 * k1R * h);
    k2R = dRdt(t + 0.5 * h, St + 0.5 * k1S * h, It + 0.5 * k1I * h, Rt + 0.5 * k1R * h);

    k3S = dSdt(t + 0.5 * h, St + 0.5 * k2S * h, It + 0.5 * k2I * h, Rt + 0.5 * k2R * h);
    k3I = dIdt(t + 0.5 * h, St + 0.5 * k2S * h, It + 0.5 * k2I * h, Rt + 0.5 * k2R * h);
    k3R = dRdt(t + 0.5 * h, St + 0.5 * k2S * h, It + 0.5 * k2I * h, Rt + 0.5 * k2R * h);

    k4S = dSdt(t + h, St + k3S * h, It + k3I, Rt + k3R);
    k4I = dIdt(t + h, St + k3S * h, It + k3I, Rt + k3R);
    k4R = dRdt(t + h, St + k3S * h, It + k3I, Rt + k3R);


    St(t + 1,1) = St(t,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t + 1,1) = It(t,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t + 1,1) = Rt(t,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

end

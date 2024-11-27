format shortG
clear
clc

%Variable Definitions
St1 = zeros(1,1);
It1 = zeros(1,1);
Rt1 = zeros(1,1);

A = 5; %Amplitude
w = 2 * pi * 365/365; %Angular Frequency
Beta0 = 0.3; %Transmission Rate IC
Gamma = 0.1; %Recovery Rate
St(1,1) = 990; %Number of susceptible individuals at time t
It(1,1) = 10; %Number of infected individuals at time t
Rt(1,1) = 0; %Number of recovered individuals at time t

N = St + It + Rt; %Total Population

h = 0.1; %Step size in days
T = 30; %Total simulation time in days (0 -> 30)

days = 1:h:T;
t1 = 1; %Counting variable

for t = 1:h:T

    Beta = Beta0 * (1 + A * sin(w * t));

    dSdt = @(t, St, It, Rt) -(Beta / N) .* St .* It;
    dIdt = @(t, St, It, Rt) (Beta / N) .* St .* It - Gamma .* It;
    dRdt = @(t, St, It, Rt) Gamma .* It;
    
    k1S = dSdt(t, St(t1,1), It(t1,1), Rt(t1,1));
    k1I = dIdt(t, St(t1,1), It(t1,1), Rt(t1,1));
    k1R = dRdt(t, St(t1,1), It(t1,1), Rt(t1,1));

    k2S = dSdt(t + 0.5 * h, St(t1,1) + 0.5 * k1S * h, It(t1,1) + 0.5 * k1I * h, Rt(t1,1) + 0.5 * k1R * h);
    k2I = dIdt(t + 0.5 * h, St(t1,1) + 0.5 * k1S * h, It(t1,1) + 0.5 * k1I * h, Rt(t1,1) + 0.5 * k1R * h);
    k2R = dRdt(t + 0.5 * h, St(t1,1) + 0.5 * k1S * h, It(t1,1) + 0.5 * k1I * h, Rt(t1,1) + 0.5 * k1R * h);

    k3S = dSdt(t + 0.5 * h, St(t1,1) + 0.5 * k2S * h, It(t1,1) + 0.5 * k2I * h, Rt(t1,1) + 0.5 * k2R * h);
    k3I = dIdt(t + 0.5 * h, St(t1,1) + 0.5 * k2S * h, It(t1,1) + 0.5 * k2I * h, Rt(t1,1) + 0.5 * k2R * h);
    k3R = dRdt(t + 0.5 * h, St(t1,1) + 0.5 * k2S * h, It(t1,1) + 0.5 * k2I * h, Rt(t1,1) + 0.5 * k2R * h);

    k4S = dSdt(t + h, St(t1,1) + k3S * h, It(t1,1) + k3I, Rt(t1,1) + k3R);
    k4I = dIdt(t + h, St(t1,1) + k3S * h, It(t1,1) + k3I, Rt(t1,1) + k3R);
    k4R = dRdt(t + h, St(t1,1) + k3S * h, It(t1,1) + k3I, Rt(t1,1) + k3R);


    St(t1 + 1,1) = St(t1,1) + (1/6) * (k1S + 2 * k2S + 2 * k3S + k4S) * h;
    It(t1 + 1,1) = It(t1,1) + (1/6) * (k1I + 2 * k2I + 2 * k3I + k4I) * h;
    Rt(t1 + 1,1) = Rt(t1,1) + (1/6) * (k1R + 2 * k2R + 2 * k3R + k4R) * h;

    St1(t1,1) = St(t1,1);
    It1(t1,1) = It(t1,1);
    Rt1(t1,1) = Rt(t1,1);

    t1 = t1 + 1;

end

figure(1)
hold on
grid on

plot(days,St1,'red','LineWidth',2)
plot(days,It1,'green','LineWidth',2)
plot(days,Rt1,'blue','linewidth',2)

title('Seasonal Influenza, ß = 0.3(1 + 5sin(wt) y = 0.1')
xlabel('Time in Days')
ylabel('People')
legend('Susceptible','Infected','Recovered')

hold off

%Q 3 has a written qustion that needs to be answered. 

fftSt = abs(fft(St1));
fftIt = abs(fft(It1));
fftRt = abs(fft(Rt1));

format shortG
clear
clc

%Variable Definitions
St1 = zeros(1,1);
It1 = zeros(1,1);
Rt1 = zeros(1,1);

A = 5; %Amplitude
w1 = 2 * pi * 365/365; %Angular Frequency
w2 = 2 * pi * 100/365; % Angular Frequency for #6
Beta0 = 0.3; %Transmission Rate IC
Gamma = 0.1; %Recovery Rate
St(1,1) = 990; %Number of susceptible individuals at time t
It(1,1) = 10; %Number of infected individuals at time t
Rt(1,1) = 0; %Number of recovered individuals at time t

N = St + It + Rt; %Total Population

h = 0.1; %Step size in days
T1 = 30; %Total simulation time in days (0 -> 30)
T2 = 30; % Total simulation time in days for #6

days = 1:h:T1;
f1 = (1/T1)*(0:length(days)-1); % frequncy vector used for plotting
days2 = (1:h:T2); 
f2 = (1/T1)*(0:length(days)-1); % frequency vector for #6
t1 = 1; %Counting variable

for t = 1:h:T1 % For loop for #5

    Beta = Beta0 * (1 + A * sin(w1 * t));

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

fftSt1 = abs(fft(St1));
fftIt1 = abs(fft(It1));
fftRt1 = abs(fft(Rt1));

m = max(fftIt1);
l = T1/2;

figure(2)
hold on
grid on

plot(f1,fftIt1,'b','LineWidth',2);
title('Fourier Analysis of Infected Cases')
legend('Discrete Fourier Transform')
ylabel('FFT')
xlabel('f')
axis([0 5, 0 m])

hold off

for t = 1:h:T1 % For loop for #6

    Beta = Beta0 * (1 + A * sin(w2 * t));

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

    St2(t1,1) = St(t1,1);
    It2(t1,1) = It(t1,1);
    Rt2(t1,1) = Rt(t1,1);

    t1 = t1 + 1;

end

fftSt2 = abs(fft(St2));
fftIt2 = abs(fft(It2));
fftRt2 = abs(fft(Rt2));

m2 = max(fftIt2);

figure(3)
hold on
grid on

plot(f2,fftIt2,'b','LineWidth',2);
title('Fourier Analysis of Infected Cases')
legend('Discrete Fourier Transform')
ylabel('FFT')
xlabel('f')
axis([0 l, 0 m2])

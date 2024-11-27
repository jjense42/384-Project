format shortG
clear
clc

%Variable Definitions
St1 = zeros(1,1);
It1 = zeros(1,1);
Rt1 = zeros(1,1);

St(1,1) = 990; %Number of susceptible individuals at time t
It(1,1) = 10; %Number of infected individuals at time t
Rt(1,1) = 0; %Number of recovered individuals at time t

N = St + It + Rt; %Total Population

h = 1; %Step size in days
T = 30; %Total simulation time in days (0 -> 100)



days = 1:h:100;

for t = 1:h:T

    Beta = 0.3; %Transmission Rate
    Gamma = 0.1; %Recovery Rate

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

    St1(t,1) = St(t,1);
    It1(t,1) = It(t,1);
    Rt1(t,1) = Rt(t,1);

end

%linear least squares with T=1:30
x=days(1:30);
y=It1';
%solve for ln I(t)=ln I(0)+kt. A1=k A0=ln I(0)
A1= (T*sum(x.*y)-(sum(x)*sum(y)))/(T*sum(x.^2)-(sum(x)^2))
A0= (sum(y)/T)-(A1*(sum(x)/T));
I0=exp(A0)
Y=[];
for t=1:T
    Y= [Y I0*exp(A1*(t))];
end

b1=(A1+Gamma)*(N/St(1,1))



%linear least squares with T=1:10
x1=days(1:10);
y1=y(1:10);

%solve for ln I(t)=ln I(0)+kt. A11=k A01=ln I(0)
A11= (10*sum(x1.*y1)-(sum(x1)*sum(y1)))/(10*sum(x1.^2)-(sum(x1)^2));
A01= (sum(y1)/10)-(A11*(sum(x1)/10));
I01=exp(A01)
Y1=[];
for t=1:10
    Y1= [Y1 I01*exp(A11*(t))];
end
b1=(A11+Gamma)*(N/St(1,1))





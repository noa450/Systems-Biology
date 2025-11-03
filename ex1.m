% Simple SIR model

%Define parameters and create value trackers
days = 600;
dt = 1;
steps = days/dt;

U = zeros(1,steps);
L = zeros(1,steps);
I = zeros(1,steps);
V = zeros(1,steps);

U(1) = 500;
L(1) = 0;
I(1) = 0;
V(1) = 10;

figure
p = 0.272  % production rate of healthy cells per day
k1 = 0.00027 % infection rate
k2 = 0.00136  % death rate of healthy and latent cells
k3 =  0.33 % death rate of infected cells
k4 = 100 % virion production per infected cell per day
k5 = 2 % virion decay rate (half day lifespan)
k6 = 0.05 % rate latent -> infectious
for step = 2:steps

    dU =  p - k2*U(step-1) - k1*U(step-1)*V(step-1);
    dL =  k1*U(step-1)*V(step-1) - k2*L(step-1) - k6*L(step-1);
    dI =  k6*L(step-1) - k3*I(step-1);
    dV =  k4*I(step-1) - k5*V(step-1);

    U(step) = U(step-1) + dU * dt;
    L(step) = L(step-1) + dL * dt;
    I(step) = I(step-1) + dI * dt;
    V(step) = V(step-1) + dV * dt;

end


hold on % This makes the plots appear on top of each other
plot(U,'b','LineWidth',2)
plot(L,'m','LineWidth',2)
plot(I,'r','LineWidth',2)
plot(V,'k','LineWidth',2)

title('HIV Model')
xlabel('Days')
ylabel('Count')
legend({'U (healthy)','L (latent)','I (infected)','V (virions)'})
xlim([0, 200])
ylim([0, 2600])
fig = gca;


% Try changing the Y axis labels to reflect the population of Israel (9.2M)


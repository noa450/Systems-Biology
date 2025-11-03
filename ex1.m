% HIV simulation - Question 3

days = 600;
dt = 1;
steps = days/dt;



% section a
U = zeros(1,steps);
L = zeros(1,steps);
I = zeros(1,steps);
V = zeros(1,steps);

U(1) = 500;
L(1) = 0;
I(1) = 0;
V(1) = 10;


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

figure
hold on % This makes the plots appear on top of each other
plot(U,'b','LineWidth',2)
plot(L,'m','LineWidth',2)
plot(I,'r','LineWidth',2)
plot(V,'k','LineWidth',2)

title('HIV Model')
xlabel('Days')
ylabel('Count/mm^3');
legend({'U (healthy)','L (latent)','I (infected)','V (virions)'})
xlim([0, 200])
ylim([0, 2600])
grid on;
fig = gca;



% section b

k6_outbreak = 0.05;% a (Latent to Infectious Activation Rate)

f_latent = 0.20;   % Fraction of infected cells that become latent
k6_dormant = 1e-10;% Effectively zero latent activation rate for Phase 1

steps = days / dt;
T_steady = 100; % Day of outbreak trigger
step_switch = T_steady / dt; 


time = 1:steps;



% Phase 1: Dormancy (t=0 to T_steady) 
for step = 1:step_switch - 1 
    
    U_curr = U(step); L_curr = L(step); I_curr = I(step); V_curr = V(step);
    
    infection_rate = k1 * U_curr * V_curr; 

    dH = p - k2 * U_curr - infection_rate;
    dL = f_latent * infection_rate - k2 * L_curr - k6_dormant * L_curr; % L accumulates (k2 = d)
    dI = (1 - f_latent) * infection_rate + k6_dormant * L_curr - k3 * I_curr; % k3 = delta
    dV = k4 * I_curr - k5 * V_curr - infection_rate; % k4=k, k5=c
    
    % Update populations for next step 
    U(step + 1) = max(0, U_curr + dH * dt);
    L(step + 1) = max(0, L_curr + dL * dt);
    I(step + 1) = max(0, I_curr + dI * dt);
    V(step + 1) = max(0, V_curr + dV * dt);
end

% Phase 2: Outbreak Triggered (t > T_steady to T_total)
for step = step_switch:steps - 1 
    
    U_curr = U(step); L_curr = L(step); I_curr = I(step); V_curr = V(step);

    infection_rate = k1 * U_curr * V_curr;
    
    dH = p - k2 * U_curr - infection_rate;
    dL = f_latent * infection_rate - k2 * L_curr - k6_outbreak * L_curr; % L is converted
    dI = (1 - f_latent) * infection_rate + k6_outbreak * L_curr - k3 * I_curr; % I spikes
    dV = k4 * I_curr - k5 * V_curr - infection_rate; 
    
    % Update populations for next step 
    U(step + 1) = max(0, U_curr + dH * dt);
    L(step + 1) = max(0, L_curr + dL * dt);
    I(step + 1) = max(0, I_curr + dI * dt);
    V(step + 1) = max(0, V_curr + dV * dt);
end

figure();
plot(time, U, 'b', 'LineWidth', 2); hold on;
plot(time, L, 'm', 'LineWidth', 2);
plot(time, I, 'r', 'LineWidth', 2);
plot(time, V, 'k', 'LineWidth', 2);

% Mark the switch point
xline(step_switch, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Outbreak Trigger');

title('HIV T-Cell Model: Outbreak After Dormancy');
xlabel('Days');
ylabel('Count/mm^3');
legend({'U (Healthy)','L (Latent)','I (Infectious)','V (Virions)'}, 'Location', 'NorthEast');
xlim([0, 200])
grid on;

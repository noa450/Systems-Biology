% HIV simulation - Question 3

days = 600;
dt = 0.01;
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
    dL =  0.2*k1*U(step-1)*V(step-1) - k2*L(step-1) - k6*L(step-1);
    dI =  0.8*k1*U(step-1)*V(step-1)+k6*L(step-1) - k3*I(step-1);
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

grid on;
fig = gca;



% section b

k6_dormant = 1e-10;% Effectively zero latent activation rate for Phase 1

T_steady = 100; % Day of outbreak trigger
step_switch = T_steady / dt; 


time = (0:steps-1)*dt;


% Phase 1: Dormancy (t=0 to T_steady) 
for step = 2:step_switch - 1 
    dU =  p - k2*U(step-1) - k1*U(step-1)*V(step-1);
    dL =  0.2*k1*U(step-1)*V(step-1) - k2*L(step-1) - k6_dormant*L(step-1);
    dI =  0.8*k1*U(step-1)*V(step-1)+k6_dormant*L(step-1) - k3*I(step-1);
    dV =  k4*I(step-1) - k5*V(step-1);

    U(step) = U(step-1) + dU * dt;
    L(step) = L(step-1) + dL * dt;
    I(step) = I(step-1) + dI * dt;
    V(step) = V(step-1) + dV * dt;
end

% Phase 2: Outbreak Triggered (t > T_steady to T_total)
for step = step_switch:steps - 1 
    
    dU =  p - k2*U(step-1) - k1*U(step-1)*V(step-1);
    dL =  0.2*k1*U(step-1)*V(step-1) - k2*L(step-1) - k6*L(step-1);
    dI =  0.8*k1*U(step-1)*V(step-1)+k6*L(step-1) - k3*I(step-1);
    dV =  k4*I(step-1) - k5*V(step-1);

    U(step) = U(step-1) + dU * dt;
    L(step) = L(step-1) + dL * dt;
    I(step) = I(step-1) + dI * dt;
    V(step) = V(step-1) + dV * dt;
end

figure();
plot(time, U, 'b', 'LineWidth', 2); hold on;
plot(time, L, 'm', 'LineWidth', 2);
plot(time, I, 'r', 'LineWidth', 2);
plot(time, V, 'k', 'LineWidth', 2);


title('HIV T-Cell Model: Outbreak After Dormancy');
xlabel('Days');
ylabel('Count/mm^3');
legend({'U (Healthy)','L (Latent)','I (Infectious)','V (Virions)'}, 'Location', 'NorthEast');
xlim([0, 200])
grid on;

%%

% question 4

epsilon_ART = 0.99; 
k4_kill_V = k4 * (1 - epsilon_ART);

k6_activated = k6 * 10;


scenarios = {'none','kill_v','activate','combined'};

Umat = zeros(length(scenarios), steps);
Lmat = zeros(length(scenarios), steps);
Imat = zeros(length(scenarios), steps);
Vmat = zeros(length(scenarios), steps);

U0 = 500; L0 = 0; I0 = 0; V0 = 10;

for s = 1:length(scenarios)
    U = zeros(1,steps); L = zeros(1,steps); I = zeros(1,steps); V = zeros(1,steps);
    U(1)=U0; L(1)=L0; I(1)=I0; V(1)=V0;
    
    switch scenarios{s}
        case 'none'
            k4_now = k4;
            k6_now = k6;
        case 'kill_v'
            k4_now = k4_kill_V;
            k6_now = k6;
        case 'activate'
            k4_now = k4;
            k6_now = k6_activated;
        case 'combined'
            k4_now = k4_kill_V;
            k6_now = k6_activated;
    end
    
    for step = 2:steps
        dU =  p - k2*U(step-1) - k1*U(step-1)*V(step-1);
        dL =  0.2*k1*U(step-1)*V(step-1) - k2*L(step-1) - k6_now*L(step-1);
        dI =  0.8*k1*U(step-1)*V(step-1)+k6_now*L(step-1) - k3*I(step-1);
        dV =  k4_now*I(step-1) - k5*V(step-1);

        U(step) = U(step-1) + dU * dt;
        L(step) = L(step-1) + dL * dt;
        I(step) = I(step-1) + dI * dt;
        V(step) = V(step-1) + dV * dt;
        
        U(step)=max(U(step),0);
        L(step)=max(L(step),0);
        I(step)=max(I(step),0);
        V(step)=max(V(step),0);
    end
    
    Umat(s,:) = U;
    Lmat(s,:) = L;
    Imat(s,:) = I;
    Vmat(s,:) = V;
end

time = (0:steps-1)*dt;
figure('Position',[100 100 900 600])
subplot(2,1,1)
hold on
plot(time, Vmat(1,:),'k','LineWidth',1.5)    % none
plot(time, Vmat(2,:),'b','LineWidth',1.5)    % ART
plot(time, Vmat(3,:),'m','LineWidth',1.5)    % activate
plot(time, Vmat(4,:),'c','LineWidth',1.5)    % combined
title('V (virions)')
xlabel('Days'); ylabel('Virions')
legend('none','ART','activate','combined')
xlim([0 200]); grid on

subplot(2,1,2)
hold on
plot(time, Lmat(1,:),'k','LineWidth',1.5)
plot(time, Lmat(2,:),'b','LineWidth',1.5)
plot(time, Lmat(3,:),'m','LineWidth',1.5)
plot(time, Lmat(4,:),'c','LineWidth',1.5)
title('L (latent)')
xlabel('Days'); ylabel('Latent cells')
legend('none','ART','activate','combined')
grid on

fprintf('\nFinal latent pool L (day %d):\n', days);
for s = 1:length(scenarios)
    fprintf('  %8s : L_final = %.4f\n', scenarios{s}, Lmat(s,end));
end

fprintf('\nRelative change in L (L_final / max(L) observed):\n');
for s = 1:length(scenarios)
    fprintf('  %8s : L_final/max(L) = %.4f\n', scenarios{s}, Lmat(s,end) / (max(Lmat(s,:)) + eps));
end
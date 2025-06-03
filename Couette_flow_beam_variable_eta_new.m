clc; close all; clear;

L_B = 64;
b = 7;
d = 5;
L = 1;
t_s = 0.1;
H = 5;
EB = 0.34; 
M = 30; 
N = 100; 


eta_min = 3;
eta_max = 10; 

eta = linspace(eta_min, eta_max, M); % Variable eta

theta0 = 0*pi/180;

h = 1e-3;

alpha_L = 0.01;   

delta_L = h/d;



n = 1.75;
C_B = 0.66;  
rho_B = 854e-9;


I_B = (1/12)*b*H^3;



for i = 1:N    
     
    Omega(i) = 1;

end


for j = 1:M

    psi_lock(j) = 0.9 / eta(j); 

    psi(j, :) = linspace(0, psi_lock(j), N); 

    psi_bar(j, :) = psi(j, :) / pi;

end


t = psi ./ Omega; 



%{
for j = 1:M

    psi_lock(j) = 0.9 / eta(j); 

end


NF = sqrt(rho_B*b*H*L_B^4/(EB*I_B));

OMEGA = 2*NF;
T = 2*pi/(1*OMEGA);


t = linspace(0,T/4,N); 


for j = 1:M
    for i = 1:N 

    psi(j,i) = psi_lock(j)*sin(OMEGA*t(i));

    Omega(j,i) = psi_lock(j)*OMEGA*cos(OMEGA*t(i));
    
    psi_bar(j,i) = psi(j,i) / pi;

    end
end

%}


for j = 1:M
    for i = 1:N
        
        psi_val = psi(j, i);
        eta_val = eta(j);
        
        % Geometric relations
        theta(j, i) = asin(eta_val * psi_val * cos(psi_val / 2)) - psi_val / 2;   

        A(j, i) = sqrt(1 - eta_val^2 * psi_val^2 * (cos(psi_val / 2))^2);   
        denomA = A(j, i) + eps;

        THETA(j, i) = (eta_val * cos(psi_val / 2) - 0.5 * eta_val * psi_val * sin(psi_val / 2)) / denomA - 0.5; 
        
        % Time derivatives
        theta_dot(j, i) = THETA(j, i) * Omega(i);

        % Derivative w.r.t psi
        term1 = (-eta_val * psi_val * cos(psi_val / 2) / 4) - eta_val * sin(psi_val / 2);
        term2 = eta_val * cos(psi_val / 2) - eta_val * sin(psi_val / 2) * psi_val / 2;
        term3 = eta_val^2 * cos(psi_val / 2) * sin(psi_val / 2) * psi_val^2 - 2 * eta_val^2 * (cos(psi_val / 2))^2 * psi_val;

        theta_dot_psi(j, i) = Omega(i) * (term1 / denomA - term2 * term3 / (2 * denomA^3 + eps));

        % R computation
        numerator = (1 / eta_val - cos(theta(j, i) - psi_val)) * sin(theta(j, i) - psi_val) * (THETA(j, i) - 1) + (psi_val / (2 * eta_val) + sin(theta(j, i) - psi_val)) * (1 / (2 * eta_val) + cos(theta(j, i) - psi_val) * (THETA(j, i) - 1));

        denominator = sqrt((1 / eta_val - cos(theta(j, i) - psi_val))^2 + (psi_val / (2 * eta_val) + sin(theta(j, i) - psi_val))^2) + eps;

        R(j, i) = numerator / denominator;

        % r_dot_bar calculation
        r_dot_bar(j, i) = R(j, i) * Omega(i) * heaviside(theta(j, i) - theta0);
        
    end
end



mu_inf_bar = 0.01;
Lambda = 0.05; 
mu_0 = 1e-6;

m = 1;


for j = 1:M
    
    epsilon(j) = h / (eta(j)*d); 
    
    a_c(j) = alpha_L * b * eta(j)*d;

end


for j = 1:M
    for i = 1:N
        
        %r_dot = r_dot_bar(j, i);
        
        mu_bar(j, i) = mu_inf_bar + (1 - mu_inf_bar) *(1 + (Lambda / epsilon(j))^2 * r_dot_bar(j,i)^2)^((m - 1) / 2);

        mu(j, i) = mu_bar(j, i) * mu_0;

        F_D(j, i) = mu(j, i) * (a_c(j) / h) * Omega(i) * eta(j)*d * R(j, i);

        F_D_bar(j, i) = mu_bar(j, i) * R(j, i) * heaviside(theta(j, i) - theta0);
        
    end
end


for j = 1:M
    for i = 1:N
        
        M_bar_linear(j, i) = psi(j, i);   

        M_bar_scales(j, i) = 12 * C_B * (L / t_s)^n * (t_s / H)^2 * (d / H) *THETA(j, i) * (theta(j, i) - theta0) * heaviside(theta(j, i) - theta0);  

        M_bar_fluidic(j, i) = 12 * Omega(i) * (mu(j, i) / EB) * (alpha_L / delta_L) * (eta(j)*d/H)^3 * R(j, i)^2 * heaviside(theta(j, i) - theta0); 
        
        M_bar(j, i) = M_bar_linear(j, i) + M_bar_scales(j, i) + M_bar_fluidic(j, i);  

    end
end


W_total = zeros(1, M);
RED = zeros(1, M);

for j = 1:M
    % Sort psi for integration
    [psi_sorted, sortIdx] = sort(psi(j, :));

    M_bar_sorted = M_bar(j, sortIdx);
    M_bar_fluidic_sorted = M_bar_fluidic(j, sortIdx);

    % Total work and RED ratio
    W_total(j) = trapz(psi_sorted, M_bar_sorted);
    RED(j) = trapz(psi_sorted, M_bar_fluidic_sorted) / W_total(j);
end


figure;
plot(eta, RED, 'o-', 'LineWidth', 2);
xlabel('\eta');
ylabel('RED = W_{fluidic} / W_{total}');
title('RED vs \eta');
grid on;

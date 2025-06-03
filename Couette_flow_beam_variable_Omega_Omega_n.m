clc; close all; clear;

% Geometry Parameters
L_B = 64;
b = 7;
l = 25; 
d = 5;
L = 1;
t_s = 0.1;
H = 5;
EB = 0.34;

eta = l / d;  

theta0 = 0*pi/180; % radians  

h = 1e-3;


alpha_L = 0.01;         
delta_L = h / d;


n = 1.75;
C_B = 0.66;  
rho_B = 854e-9;


N = 200; 
numOmega = 30;


I_B = (1/12) * b * H^3;
psi_lock = 0.9 / eta;


%NF = sqrt(rho_B * b * H * L_B^4 / (EB * I_B));
NF = pi^2*sqrt((EB*I_B)/(rho_B*b*H*L_B^4));

ratio = linspace(0.5, 2, numOmega);


OMEGA = zeros(1, numOmega);
T = zeros(1, numOmega);
t = zeros(numOmega, N); 


for j = 1:numOmega

    OMEGA(j) = ratio(j) * NF;
    T(j) = 2 * pi / OMEGA(j);
    t(j, :) = linspace(0, T(j)/4, N); 

end


psi = zeros(numOmega, N);
Omega = zeros(numOmega, N);
psi_bar = zeros(numOmega, N);

for j = 1:numOmega
    psi(j, :) = psi_lock * sin(OMEGA(j) * t(j, :));
    Omega(j, :) = psi_lock * OMEGA(j) * cos(OMEGA(j) * t(j, :));
    psi_bar(j, :) = psi(j, :) / pi; 
end


theta = zeros(numOmega, N);
A = zeros(numOmega, N);
THETA = zeros(numOmega, N);
theta_dot = zeros(numOmega, N);
theta_dot_psi = zeros(numOmega, N);
R = zeros(numOmega, N);
r_dot_bar = zeros(numOmega, N);

for j = 1:numOmega
    for i = 1:N
        psi_val = psi(j, i);
        
        theta(j, i) = asin(eta * psi_val * cos(psi_val / 2)) - psi_val / 2;   

        A(j, i) = sqrt(1 - eta^2 * psi_val^2 * (cos(psi_val / 2))^2);   

        denomA = A(j, i) + eps;
        
        THETA(j, i) = (eta * cos(psi_val / 2) - 0.5 * eta * psi_val * sin(psi_val / 2)) / denomA - 0.5; 

        theta_dot(j, i) = THETA(j, i) * Omega(j, i);

        % theta_dot_psi derivative calculation
        term1 = (-eta * psi_val * cos(psi_val / 2) / 4) - eta * sin(psi_val / 2);
        term2 = eta * cos(psi_val / 2) - eta * sin(psi_val / 2) * psi_val / 2;
        term3 = eta^2 * cos(psi_val / 2) * sin(psi_val / 2) * psi_val^2 ...
                - 2 * eta^2 * (cos(psi_val / 2))^2 * psi_val;

        theta_dot_psi(j, i) = Omega(j, i) * ...
            (term1 / denomA - term2 * term3 / (2 * denomA^3 + eps));

        % R calculation with numerator and denominator
        numerator = (1 / eta - cos(theta(j, i) - psi_val)) * sin(theta(j, i) - psi_val) * (THETA(j, i) - 1) + (psi_val / (2 * eta) + sin(theta(j, i) - psi_val)) * (1 / (2 * eta) + cos(theta(j, i) - psi_val) * (THETA(j, i) - 1));

        denominator = sqrt((1 / eta - cos(theta(j, i) - psi_val))^2 + (psi_val / (2 * eta) + sin(theta(j, i) - psi_val))^2) + eps;

        R(j, i) = numerator / denominator;

        % r_dot_bar calculation
        r_dot_bar(j, i) = R(j, i) * Omega(j, i) * heaviside(theta(j, i) - theta0);
    end
end



mu_inf_bar = 0.01;
Lambda = 0.05; 
mu_0 = 1e-6;
epsilon = h /l; 
m = 1; 


mu_bar = zeros(numOmega, N);
mu = zeros(numOmega, N);
F_D_bar = zeros(numOmega, N);

for j = 1:numOmega
    for i = 1:N
        r_dot = r_dot_bar(j, i);

        mu_bar(j, i) = mu_inf_bar + (1 - mu_inf_bar) * ...
            (1 + (Lambda / epsilon)^2 * r_dot^2)^((m - 1) / 2);

        mu(j, i) = mu_bar(j, i) * mu_0;

        F_D_bar(j, i) = mu_bar(j, i) * R(j, i) * heaviside(theta(j, i) - theta0);
    end
end


M_bar_linear = zeros(numOmega, N);
M_bar_scales = zeros(numOmega, N);
M_bar_fluidic = zeros(numOmega, N);
M_bar = zeros(numOmega, N);

for j = 1:numOmega
    for i = 1:N
        M_bar_linear(j, i) = psi(j, i);

        M_bar_scales(j, i) = 12 * C_B * (L / t_s)^n * (t_s / H)^2 * (d / H) * THETA(j, i) * (theta(j, i) - theta0) * heaviside(theta(j, i) - theta0);

        M_bar_fluidic(j, i) = 12 * Omega(j, i) * (mu(j, i) / EB) * (alpha_L / delta_L) * (l / H)^3 * R(j, i)^2 * heaviside(theta(j, i) - theta0);

        M_bar(j, i) = M_bar_linear(j, i) + M_bar_scales(j, i) + M_bar_fluidic(j, i);

    end
end


W_total = zeros(1, numOmega);
RED = zeros(1, numOmega);

for j = 1:numOmega
    
    [psi_sorted, sortIdx] = sort(psi(j, :));
    
    M_bar_sorted = M_bar(j, sortIdx);
    M_bar_fluidic_sorted = M_bar_fluidic(j, sortIdx);
    
    W_total(j) = trapz(psi_sorted, M_bar_sorted);
    
    RED(j) = trapz(psi_sorted, M_bar_fluidic_sorted) / W_total(j);

end
% Updated script: variable alpha_L with slow‑then‑fast spacing
clc; clear; close all;

%% --- GEOMETRY & CONSTANTS ---
b    = 7;       % mm
l    = 25;      % mm
d    = 5;       % mm
L    = 1;       % mm
t_s  = 0.1;     % mm
H    = 5;       % mm
EB   = 0.34;    % normalized modulus

% single film thickness
h       = 1e-3;          % mm
delta_L = h / d;         % dimensionless gap
epsilon = h / l;         % dimensionless film thickness

%% --- NONLINEAR α_L SPACING ---
alpha_min = 0.001;
alpha_max = 0.1;
N_alpha   = 30;
p         = 3;           % exponent >1 → slow at start, fast at end

u         = linspace(0,1,N_alpha);
alpha_L   = alpha_min + (alpha_max - alpha_min) * (u.^p);
J         = numel(alpha_L);
a_c       = alpha_L * b * l;  % contact area factor

%% --- FLUID PARAMETERS ---
mu_inf_bar = 0.01;
Lambda     = 0.05;
mu_0       = 1e-6;   % Pa·s
m          = 1;

%% --- KINEMATIC PRECOMPUTATIONS ---
eta      = l / d;
theta0   = 58*pi/180;                % rad
psi_lock = 0.9 / eta;        % engagement limit
N        = 100;
psi      = linspace(0, psi_lock, N);
Omega    = ones(size(psi));  % rad/s

% scale rotation kinematics
theta   = asin(eta .* psi .* cos(psi/2)) - psi/2;
A       = sqrt(1 - (eta .* psi).^2 .* cos(psi/2).^2);
THETA   = (eta .* cos(psi/2) - (eta.*psi/2).*sin(psi/2)) ./ A - 0.5;

% contact radius rate
num1      = (1/eta - cos(theta-psi)).*sin(theta-psi).*(THETA-1) ...
           + (psi/(2*eta) + sin(theta-psi)) ...
           .* (1/(2*eta) + cos(theta-psi).*(THETA-1));
den1      = sqrt((1/eta - cos(theta-psi)).^2 + (psi/(2*eta)+sin(theta-psi)).^2);
R         = num1 ./ den1;
r_dot_bar = R .* Omega .* heaviside(theta - theta0);

%% --- MOMENT COMPONENTS (ψ‑dependent only) ---
M_bar_linear  = psi;
M_bar_scales  = 12 * 0.66 * (L/t_s)^1.75 * (t_s/H)^2 * (d/H) ...
                .* THETA .* (theta - theta0) .* heaviside(theta - theta0);
M_bar_elastic = M_bar_linear + M_bar_scales;

%% --- PREALLOCATE MATRICES ---
mu_bar        = zeros(J, N);
mu            = zeros(J, N);
F_D           = zeros(J, N);
F_D_bar       = zeros(J, N);
M_bar_fluidic = zeros(J, N);
M_total       = zeros(J, N);
W_total       = zeros(1, J);
W_fluid       = zeros(1, J);

%% --- LOOP OVER α_L VALUES ---
for j = 1:J
    % shear‐rate dependent viscosity
    mu_bar(j,:) = mu_inf_bar + (1 - mu_inf_bar) ...
                  .* (1 + (Lambda/epsilon)^2 .* r_dot_bar.^2).^((m-1)/2);
    mu(j,:) = mu_bar(j,:) * mu_0;
    
    % drag force (optional)
    F_D(j,:)     = mu(j,:) .* (a_c(j)/h) .* Omega .* l .* R;
    F_D_bar(j,:) = mu_bar(j,:) .* R .* heaviside(theta - theta0);
    
    % fluidic bending moment
    M_bar_fluidic(j,:) = 12 .* Omega ...
                         .* (mu(j,:)/EB) ...
                         * (alpha_L(j)/delta_L) ...
                         * (l/H)^3 ...
                         .* R.^2 ...
                         .* heaviside(theta - theta0);
    
    % total moment per ψ
    M_total(j,:) = M_bar_elastic + M_bar_fluidic(j,:);
    
    % integrate work
    W_total(j) = trapz(psi, M_total(j,:));
    W_fluid(j) = trapz(psi, M_bar_fluidic(j,:));
end

% energy‐dissipation ratio vs alpha_L
RED = W_fluid ./ W_total;

%% --- PLOT RED vs α_L (log–log) ---
figure;
loglog(alpha_L, RED, '.', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('\alpha_L');
ylabel('RED = W_{fluid}/W_{total}');
title('RED vs \alpha_L (slow→fast sampling)');
grid on;

% Updated script: variable delta_L
clc; clear; close all;


L_B    = 64;        % mm (unused in this snippet)
b      = 7;         % mm
l      = 25;        % mm
d      = 5;         % mm
L      = 1;         % mm
t_s    = 0.1;       % mm
H      = 5;         % mm
EB     = 0.34;      % normalized


%h      = linspace(1e-4,1e-2,50);   % mm

h_min = 1e-4;       % starting film thickness (mm)
h_max = 1e-2;       % ending film thickness (mm)
N_h   = 30;         % total number of points
p     = 3;          % exponent >1 → slow at start, fast at end

u = linspace(0, 1, N_h);
h = h_min + (h_max - h_min) * (u.^p);


delta_L = h / d;                    % dimensionless gap parameter
epsilon = h / l;                    % dimensionless film thickness


alpha_L = 0.01;                     % lubrication area ratio


n         = 1.75;
C_B       = 0.66;
rho_B     = 854e-9;   % assume kg/mm^3 or normalized
mu_inf_bar= 0.01;
Lambda    = 0.05;
mu_0      = 1e-6;     % Pa·s
m         = 1;


eta      = l / d;
theta0   = 0*pi/180;                     % rad
psi_lock = 0.9 / eta;            % engagement limit
N        = 100;
psi      = linspace(0, psi_lock, N);
Omega    = 1*ones(size(psi));      % constant angular rate



theta   = asin(eta .* psi .* cos(psi/2)) - psi/2;
A       = sqrt(1 - (eta .* psi).^2 .* cos(psi/2).^2);
THETA   = (eta .* cos(psi/2) - (eta.*psi/2).*sin(psi/2)) ./ A - 0.5;


num1    = (1/eta - cos(theta-psi)).*sin(theta-psi).*(THETA-1) ...
         + (psi/(2*eta) + sin(theta-psi)).*(1/(2*eta) + cos(theta-psi).*(THETA-1));
den1    = sqrt((1/eta - cos(theta-psi)).^2 + (psi/(2*eta)+sin(theta-psi)).^2);
R       = num1 ./ den1;
r_dot_bar = R .* Omega .* heaviside(theta - theta0);


M_bar_linear = psi;
M_bar_scales = 12 * C_B * (L/t_s)^n * (t_s/H)^2 * (d/H) .* THETA .* (theta - theta0) .* heaviside(theta - theta0);


J = numel(delta_L);
mu_bar        = zeros(J, N);
mu            = zeros(J, N);
M_bar_fluidic = zeros(J, N);
M_total       = zeros(J, N);
W_total       = zeros(1, J);
W_elastic     = zeros(1, J);
W_fluid       = zeros(1, J);


for j = 1:J

    mu_bar(j,:) = mu_inf_bar + (1 - mu_inf_bar)*(1 + ( (Lambda/epsilon(j)).^2 ) .* r_dot_bar.^2 ).^((m-1)/2);
    mu(j,:)     = mu_bar(j,:) * mu_0;
    

    M_bar_fluidic(j,:) = 12 .* Omega .* (mu(j,:)/EB) * (alpha_L / delta_L(j)) * (l/H)^3 .* R.^2 .* heaviside(theta - theta0);
                     

    M_total(j,:) = M_bar_linear + M_bar_scales + M_bar_fluidic(j,:);
    

    W_total(j)   = trapz(psi, M_total(j,:));
    W_elastic(j) = trapz(psi, M_bar_linear + M_bar_scales);
    W_fluid(j)   = trapz(psi, M_bar_fluidic(j,:));

end


RED = W_fluid ./ W_total;


X1 = log(delta_L);
X2 = log(alpha_L);
X3 = log(RED);


%% --- PLOT RED vs δL ---
figure;
loglog(delta_L, RED, '.');
%plot(delta_L, RED, 'LineWidth', 2);
xlabel('\delta_L = h/d');
ylabel('RED = W_{fluid}/W_{total}');
title('Regime map: RED vs \delta_L');
grid on;


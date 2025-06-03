%% --- SETUP & GEOMETRY ---
clc; clear; close all;

% Geometry / substrate properties
L_B    = 64;        % mm (unused here)
b      = 7;         % mm
l      = 25;        % mm
d      = 5;         % mm
L      = 1;         % mm
t_s    = 0.1;       % mm
H      = 5;         % mm
EB     = 0.34;      % normalized modulus

%% --- FILM THICKNESS / δ_L VECTOR ---
h_min = 1e-4;       
h_max = 1e-2;       
N_h   = 30;         
p_h   = 3;          

u_h   = linspace(0,1,N_h);
h     = h_min + (h_max - h_min)*(u_h.^p_h);

delta_L   = h / d;           % (1×N_h)
epsilon   = h / l;           % (1×N_h)
N_delta   = numel(delta_L);

%% --- α_L VECTOR ---
alpha_min = 0.001;
alpha_max = 0.1;
N_alpha   = 30;
p_a       = 3;           

u_a       = linspace(0,1,N_alpha);
alpha_L   = alpha_min + (alpha_max - alpha_min)*(u_a.^p_a);
N_alpha   = numel(alpha_L);

a_c = alpha_L * b * l;      % (1×N_alpha)

%% --- FLUID & SCALE PARAMETERS ---
n           = 1.75;
C_B         = 0.66;
rho_B       = 854e-9;       % kg/mm^3 (or normalized)
mu_inf_bar  = 0.01;
Lambda      = 0.05;
mu_0        = 1e-6;         % Pa·s
m           = 1;

eta         = l / d;
theta0      = 0;            % initial inclination
psi_lock    = 0.9/eta;      
N_psi       = 100;
psi         = linspace(0,psi_lock,N_psi);
Omega       = ones(size(psi));  % constant angular rate

theta   = asin(eta.*psi.*cos(psi/2)) - psi/2;
A       = sqrt(1 - (eta.*psi).^2 .* cos(psi/2).^2);
THETA   = (eta.*cos(psi/2) - (eta.*psi/2).*sin(psi/2))./A - 0.5;

% scale‐only moment
M_bar_linear = psi;
M_bar_scales = 12*C_B*(L/t_s)^n*(t_s/H)^2*(d/H).*THETA.*(theta-theta0).*heaviside(theta-theta0);

% rotation rate term
num1    = (1/eta - cos(theta-psi)).*sin(theta-psi).*(THETA-1) ...
         + (psi/(2*eta) + sin(theta-psi)).*(1/(2*eta) + cos(theta-psi).*(THETA-1));
den1    = sqrt((1/eta - cos(theta-psi)).^2 + (psi/(2*eta)+sin(theta-psi)).^2);
R       = num1./den1;
r_dot_bar = R .* Omega .* heaviside(theta - theta0);

%% --- PREALLOCATE MATRICES ---
% Dimensions: (α-index, δ-index, psi-index)
mu_bar        = zeros(N_alpha, N_delta, N_psi);
mu            = zeros(N_alpha, N_delta, N_psi);
M_bar_fluidic = zeros(N_alpha, N_delta, N_psi);
M_total       = zeros(N_alpha, N_delta, N_psi);

W_total       = zeros(N_alpha, N_delta);
W_elastic     = zeros(N_alpha, N_delta);
W_fluid       = zeros(N_alpha, N_delta);

%% --- DOUBLE LOOP OVER α_L AND δ_L ---
for i = 1:N_alpha
    for j = 1:N_delta
        % local scalars
        eps_ij = epsilon(j);
        dL_ij  = delta_L(j);
        aL_ij  = alpha_L(i);
        
        % viscosity
        mu_bar(i,j,:) = mu_inf_bar + (1-mu_inf_bar) * (1 + ((Lambda/eps_ij)^2).*r_dot_bar.^2).^((m-1)/2);
        mu(i,j,:)     = squeeze(mu_bar(i,j,:)) * mu_0;
        
        % fluidic moment
        %M_bar_fluidic(i,j,:) = 12 .* Omega .* (squeeze(mu(i,j,:))/EB) * (aL_ij/dL_ij) * (l/H)^3 .* R.^2 .* heaviside(theta - theta0);

        M_bar_fluidic(i,j,:) = 12 * (aL_ij/delta_L(j)) * (l/H)^3 .* ( Omega .* (squeeze(mu(i,j,:))'/EB) .* R.^2 .* heaviside(theta - theta0) );
                   
        % total moment
        M_total(i,j,:) = M_bar_linear + M_bar_scales + squeeze(M_bar_fluidic(i,j,:))';
        
        % work integrals
        W_total(i,j)   = trapz(psi, M_total(i,j,:));
        W_elastic(i,j) = trapz(psi, M_bar_linear + M_bar_scales);
        W_fluid(i,j)   = trapz(psi, squeeze(M_bar_fluidic(i,j,:)));
    end
end

% Relative energy dissipation
RED = W_fluid ./ W_total;

X1 = log(alpha_L);
X2 = log(delta_L);




%% --- Surface heat‐map of RED on log–log axes ---
figure('Name','RED vs \log(\alpha_L), \log(\delta_L)','Color','w');
colormap('jet');

% build grid of log‐10 values
[DL_grid, AL_grid] = meshgrid(log10(delta_L), log10(alpha_L));

% plot
surf(DL_grid, AL_grid, RED, 'EdgeColor','none');
shading interp;
view(2);                    % top–down view

% styling
axesHandle = gca;
set(axesHandle, ...
    'FontSize',18, ...
    'FontName','Times New Roman', ...
    'Box','on');

% labels & title
xlabel('$\log(\delta_L)$', ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'FontName','Times New Roman');
ylabel('$\log(\alpha_L)$', ...
    'Interpreter','latex', ...
    'FontSize',18, ...
    'FontName','Times New Roman', ...
    'FontWeight','bold');
title('RED', ...
    'Interpreter','latex', ...
    'FontSize',14, ...
    'FontName','Times New Roman');

% axis limits
xlim([min(log10(delta_L)), max(log10(delta_L))]);
ylim([min(log10(alpha_L)), max(log10(alpha_L))]);

% DL_grid and AL_grid are meshgrids of log10(delta_L), log10(alpha_L):
%Z_iso = DL_grid + AL_grid;    

% then draw contours Z_iso = C_i:
%iso_levels = linspace(min(Z_iso(:)), max(Z_iso(:)),5);
%contour(DL_grid, AL_grid, Z_iso, iso_levels, '--k');
%clabel(contour(...), 'FontSize',12);
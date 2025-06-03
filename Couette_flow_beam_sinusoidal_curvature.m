%fixed eta
clc; close all; clear all;

L_B = 64;
b = 7;

l = 25; 
d = 5;
L = 1;
t_s = 0.1;
H = 5;
EB = 0.34; 
 

%%%------- GEOMETRY BLOCK ------ 
eta = l/d;        
theta0 = 0*pi/180;     

h = 1e-3;

alpha_L = 0.01;      
delta_L = h/d;  


%%%------- MATERIAL PARAMETER BLOCK ---------
    
n = 1.75;C_B = 0.66;  

rho_B = 854e-9;

N = 200; 


I_B = (1/12)*b*H^3;


psi_lock = 0.9/eta;


%NF1 = sqrt(rho_B*b*H*L_B^4/(EB*I_B));


NF = pi^2*sqrt((EB*I_B)/(rho_B*b*H*L_B^4));


OMEGA = 1*NF;
T = 2*pi/(1*OMEGA);


t = linspace(0,T,N); 


for i = 1:N 

    psi(i) = psi_lock*sin(OMEGA*t(i));

    Omega(i) = psi_lock*OMEGA*cos(OMEGA*t(i));
     

end

psi_bar = psi/3.1416; 





for i = 1:N
    
    theta(i) = asin(eta*psi(i)*cos(psi(i)/2))-psi(i)/2;   

    A(i) = sqrt(1-eta^2*psi(i)^2*(cos(psi(i)/2)^2));   
    
    THETA(i) = (eta*cos(psi(i)/2)-1/2*eta*psi(i)*sin(psi(i)/2))/A(i)-1/2; 
    
    %R1(i) =   sin(psi(i))*((0.5+THETA(i))/(eta^2*psi(i)^2*cos(psi(i)/2)^2)-1) - cos(psi(i))*A(i)/(eta*psi(i)*cos(psi(i)/2));  

    theta_dot(i) = THETA(i)*Omega(i);

    theta_dot_psi(i) = Omega(i)*(((-eta*psi(i)*cos(psi(i)/2)/4)-eta*sin(psi(i)/2))/A(i)-((eta*cos(psi(i)/2)- eta*sin(psi(i)/2)*psi(i)/2)*((eta^2*cos(psi(i)/2)*sin(psi(i)/2)*psi(i)^2-2*eta^2*(cos(psi(i)/2))^2*psi(i)))/(2*A(i)^3)));
    
    numerator1(i) = (1/eta-cos(theta(i)-psi(i)))*sin(theta(i)-psi(i))*(THETA(i)-1)+(psi(i)/(2*eta)+sin(theta(i)-psi(i)))*(1/(2*eta)+cos(theta(i)-psi(i))*(THETA(i)-1));
    
    denominator1(i) = sqrt((1/eta-cos(theta(i)-psi(i)))^2+(psi(i)/(2*eta)+sin(theta(i)-psi(i)))^2);
    
    R(i) = numerator1(i) / denominator1(i);

    r_dot_bar(i) = R(i)*Omega(i)*heaviside(theta(i)-theta0);
    
end


mu_inf_bar = 0.01;
Lambda = 0.05; 
mu_0 = 1e-6;
epsilon = h/l; 
m = 1;


for i = 1:N

    %mu_bar(i) = (mu_inf_bar + (1-mu_inf_bar)*(1 + (Lambda^2/(h/l)^2)*(R(i)^2)*(Omega(i)^2))^((m-1)/2))*heaviside(theta(i)-theta0);
     
    mu_bar(i) = mu_inf_bar + (1-mu_inf_bar)*(1 + (Lambda/epsilon)^2*(r_dot_bar(i))^2)^((m-1)/2);
    
    X1(i) = ((r_dot_bar(i))^2)^((m-1)/2); 

    mu(i) = mu_bar(i)*mu_0;
    
    F_D_bar(i) = mu_bar(i).*R(i)*heaviside(theta(i)-theta0); 

end


for i = 1:N

    M_bar_linear(i) = psi(i);   
    
    M_bar_scales(i) = 12*C_B*(L/t_s)^n*(t_s/H)^2*(d/H)*THETA(i)*(theta(i)-theta0)*heaviside(theta(i)-theta0);  
    
    M_bar_fluidic(i) = 12*Omega(i)*(mu(i)/EB)*alpha_L/delta_L*(l/H)^3*R(i)^2*heaviside(theta(i)-theta0); 
    
    M_bar(i) = M_bar_linear(i) + M_bar_scales(i) + M_bar_fluidic(i);  
    
end


Dissipative_energy = trapz(psi(1:N), M_bar_fluidic(1:N));

W_total = trapz(psi(1:N), M_bar(1:N));

RED = Dissipative_energy / W_total;





figure(1); 
grid on; hold on;
plot(t, M_bar, 'k.-'); 
xlabel('t (sec)', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$\bar{M}$', 'interpreter', 'latex', 'FontSize', 16);
title('Moment-Curvature for $\eta = 6$, $\theta_0 = 0$');



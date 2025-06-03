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


alpha_L = 0.01; %0.5/l;         
delta_L = h/d;  


%%%------- MATERIAL PARAMETER BLOCK ---------
n = 1.75;C_B = 0.66;  

rho_B = 854e-9;

N = 100; 


I_B = (1/12)*b*H^3;


psi_lock = 0.9/eta;


psi = linspace(0,psi_lock,N); 


psi_bar = psi/3.1416;


for i = 1:N    
     
    Omega(i) = 1;

end



t = psi./Omega; 



for i = 1:N
    
    theta(i) = asin(eta*psi(i)*cos(psi(i)/2))-psi(i)/2;   

    A(i) = sqrt(1-eta^2*psi(i)^2*(cos(psi(i)/2)^2));   
    
    THETA(i) = (eta*cos(psi(i)/2)-1/2*eta*psi(i)*sin(psi(i)/2))/A(i)-1/2; 
    
    %R(i) =   sin(psi(i))*((0.5+THETA(i))/(eta^2*psi(i)^2*cos(psi(i)/2)^2)-1) - cos(psi(i))*A(i)/(eta*psi(i)*cos(psi(i)/2));  

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

a_c = alpha_L*b*l;

K_theta = (b)*EB*0.66*t_s^2*(L/t_s)^1.75;


for i = 1:N

    %mu_bar1(i) = mu_inf_bar + (1-(mu_inf_bar))*(1 +(Lambda^2/(h/l)^2)*(R(i)^2)*(Omega(i)^2))^((m-1)/2);
        
    mu_bar(i) = mu_inf_bar + (1-mu_inf_bar)*(1 + (Lambda/epsilon)^2*(r_dot_bar(i))^2)^((m-1)/2);

    mu(i) = mu_bar(i)*mu_0;
    
    F_D(i) = mu(i)*(a_c/h)*Omega(i)*l*R(i); 

    %F_D_bar(i) = mu_bar(i).*R(i)*heaviside(theta(i)-theta0); 

    F_D_bar(i) = F_D(i)*d/(K_theta);
    %F_D_bar1(i) = F_D(i)*l/(K_theta*(theta(i)-theta0));

    Wi(i) = Lambda*l*r_dot_bar(i)/epsilon;

end



for i = 1:N

    M_bar_linear(i) = psi(i);   
    
    M_bar_scales(i) = 12*C_B*(L/t_s)^n*(t_s/H)^2*(d/H)*THETA(i)*(theta(i)-theta0)*heaviside(theta(i)-theta0);  
    
    M_bar_fluidic(i) = 12*Omega(i)*(mu(i)/EB)*(alpha_L/delta_L)*(l/H)^3*R(i)^2*heaviside(theta(i)-theta0); 
    
    M_bar(i) = M_bar_linear(i) + M_bar_scales(i) + M_bar_fluidic(i);  
     
    M_bar_elastic(i) = M_bar_linear(i) + M_bar_scales(i);
end



W_total = trapz(psi(1:end), M_bar(1:end));
RED = trapz(psi(1:end), M_bar_fluidic(1:end)) / W_total;

W_elastic = trapz(psi(1:end), M_bar_elastic(1:end));
W_fluid = trapz(psi(1:end), M_bar_fluidic(1:end));





figure(1); 
grid on; hold on;
plot(t, psi_bar, 'k.-'); 
xlabel('t (sec)', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$\bar{\psi}$', 'interpreter', 'latex', 'FontSize', 16);

figure(2); 
grid on; hold on;
plot(psi_bar, Wi, 'k.-'); 
xlabel('$\bar{\psi}$', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$Wi$', 'interpreter', 'latex', 'FontSize', 16);

figure(3); 
grid on; hold on;
plot(r_dot_bar, mu_bar, 'k.-'); 
xlabel('$\dot{\bar{r}}$', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$\bar{\mu}$', 'interpreter', 'latex', 'FontSize', 16);


figure(4); 
grid on; hold on;
plot(psi_bar, F_D_bar, 'k.-'); 
xlabel('$\bar{\psi}$', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$\pi_D$', 'interpreter', 'latex', 'FontSize', 16);


figure(5); 
grid on; hold on;
plot(psi_bar, M_bar, 'k.-'); 
xlabel('$\bar{\psi}$', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$\bar{M}$', 'interpreter', 'latex', 'FontSize', 16);
title('Moment-Curvature for $\eta = 6$, $\theta_0 = 0$');

figure(6); grid on; hold on;
plot(psi_bar, M_bar_fluidic, 'k.-'); 
xlabel('$\bar{\psi}$', 'interpreter', 'latex', 'FontSize', 16); 
ylabel('$\bar{M}$', 'interpreter', 'latex', 'FontSize', 16);
title('Fluidic Damping for $\eta = 6$, $\theta_0 = 0$');



    %C1(i) = Omega(i)*mu(i);

    %C2(i) = Omega(i)*mu(i)*R(i)^2;

    %M_bar_fluidic1(i) = (l*d)/(EB*I_B)*F_D(i)*R(i); 


%figure(6);
%subplot(2,2,3)
%grid minor;hold on;
%plot(Omega(2:N),B_bar(2:N),'k.-');
%plot(psi(2:N/2),B_bar(2:N/2)/1e3,'k.-');
%plot(psi((N/2):end-1),B_bar((N/2):end-1)/1e3,'r.-');
%xlabel('$\Omega$ $(rad)$','interpreter','latex','FontSize',16)
%ylabel('\bar{B}')
%title('Bending rigidity-Curvature for \eta = 6, \theta_0=0')


%{
rho_b = linspace(5000,20000,16);
eta1 = linspace(3,10,500);

AA = xlsread('Book2.xlsx');

figure(7)

[X,Y] = meshgrid(eta1,rho_b);

surf(X',Y',AA)
shading interp 
colorbar
xlabel('\eta')
ylabel('\rho')
xlim([3 10])
%}


%{
for i = 1:N

    M_bar_linear(i) = eta*psi(i);   
    
    M_bar_scales(i) = 12*C_B*A_LD^(n-2)*A_LH^3*A_lL*THETA(i)*(theta(i)-theta0)*heaviside(theta(i)-theta0);  
    
    M_bar_micro_inerta(i) = (l/(E*I))*I1*theta_dot(i)*theta_dot_psi(i);
    
    M_bar_fluidic(i) = 12*eta*Omega(i)*tau*alpha_L/delta_L*A_LH^3*A_lL^3*R(i)^2*heaviside(theta(i)-theta0); 
    
    M_bar(i) = M_bar_linear(i) + M_bar_scales(i) + M_bar_fluidic(i)+M_bar_micro_inerta(i);  
    
end
%}

%{
for i = 1:N

    Mz_bar_linear(i) =(L_B/d)*psi(i);   
    
    Mz_bar_scales(i) = 12*C_B*(L/t_s)^n*(t_s/H)^2*(L_B/H)*(beta/(beta))*THETA(i)*(theta(i)-theta0)*heaviside(theta(i)-theta0);  
    
    Mz_bar_fluidic(i) = 12*(L_B/d)*(beta/(beta))*Omega(i)*tau*(alpha_L/delta_L)*(l/H)^3*R(i)^2*heaviside(theta(i)-theta0); 
    
    Mz_bar(i) = Mz_bar_linear(i) + Mz_bar_scales(i) + Mz_bar_fluidic(i);  
    
end
%}
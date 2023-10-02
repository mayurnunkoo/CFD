%% General Information
%Project: Finite Volume Method, CDS/upwind
%Name: Mayur Nunkoo
%Education: Aerospace Engineering
%Capstone: BEFAV
%Date: 2023-10-02
%Revision Number: 0
%% Clean
clc
clear
close all

%% General Parameters Definition
%length of domain
L = input("Enter length of domain: ");
%number of internal grid points
n = input("Enter the number of grid points: ");
%velocity
u = 1.5;
%centroid - centroid distance
delx = L/n;
%Diffusivity
gamma = 0.1;
%density
rho = 1;
%% Interior Parameters Definition
De = gamma/delx;
Dw = De;
Fe = rho*u;
Fw = Fe;
ae = De - 0.5*Fe
aw = Dw + 0.5*Fw
ap = ae + aw;
%% Boundary Coefficients
%Boundary conditions
phi_left = 1;
phi_right = 0;
%left coefficients
aeb_left = 0.5*Fe - De;
awb_left = -Fw - (gamma/(delx/2));
awb_left_phi = awb_left*phi_left;
apb = 0.5*Fe + De + (gamma/(delx/2));
%right coefficients
aeb_right = Fe - 2*De;
awb_right = -Fw*0.5 - Dw;
aeb_right_phi = aeb_right*phi_right;
apb_right =-Fw*0.5 + 2*De + Dw;
%% Matrix creation
%coefficient matrix
A = zeros(n);
for i = 2:n-1
    for j = i:n
       if j == i
            A(i,j-1) = -aw
            A(i,j) = ap
            A(i,j+1) = -ae
       else
            break
       end
    end
    
end
%Populating boundary coefficients
A(1,1) = apb
A(1,2) = aeb_left
A(end,end-1) = awb_right
A(end,end)= apb_right
%% Solution matrix
for i = 1:n
    if i == 1
        sol(i) = -awb_left_phi
    else
        sol(i) = 0
        
    end
end
%Transpose
solt = sol'
%% Numerical Solution
sol_computational = inv(A)*solt
%% Analytical Solution
phi0 = phi_left;
phiL = phi_right;
A = (phiL-phi0);
d = (rho*u*L)/gamma;
%% x domain creation
xdomain (1) = 1/(2*n);
for i = 1:n-1
    xdomain(i+1) = (1/(2*n)) + (1/n)*i
end

for i = 1:n
    c(i) = (rho*u*xdomain(i))/gamma;
    sol_analytical(i) = A*((exp(c(i))-1)/(exp(d)-1))+phi_left;
end
%% Plotting
figure(1)
plot(xdomain, sol_analytical,'--*')
hold on
plot(xdomain, sol_computational, '-o')
xlabel('x(m)')
ylabel('\phi')
legend('analytical','computational(CDS)')
title("Analytical vs Computational Solution (High Peclet Number)")

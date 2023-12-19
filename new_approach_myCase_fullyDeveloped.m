clc;
clear;

W = 100e-6;
nions = 4;
Ny = 51;

y = linspace(0,W,Ny).';
dy = y(2) - y(1);
y_Jy = linspace(0-dy/2,W+dy/2,Ny+1).';

D = [1.330e-9 2.030e-9 9.310e-9 1.180e-9]; % Na, Cl, H, HCO3
Z = [1 -1 1 -1];

dt = dy^2 / max(D) / 5;
tfinal = 30;
tol = 1e-12;

rho_NaCl_in = 1e-3;
rho_H_high = 1;
rho_H_low = 1e-3;

rho = zeros(Ny,nions);
phi = zeros(Ny,1);
temp = zeros(Ny,1);
Jy = zeros(Ny+1,nions);

% Initialize the rho fields
rho(:,1) = rho(:,1) + rho_NaCl_in;
rho(:,2) = rho(:,1);
rho(:,3) = (rho_H_high-rho_H_low)/W*y + rho_H_low;
rho(:,4) = rho(:,3);

dphi = phi;
rho_new = rho;

counter2 = 1;
for t = 0:dt:tfinal
    
    g = zeros(Ny,1);
    sigma = zeros(Ny,1);
    sumrhoz = zeros(Ny,1);
    for k = 1:nions
        g = g + Z(k) * D(k) * rho(:,k);
        sigma = sigma + Z(k)^2 * D(k) * rho(:,k);
        sumrhoz = sumrhoz + Z(k) * rho(:,k);
    end
    
    % Find the electric potential by relaxation
    counter = 1;
    while true
        for j = 1:Ny
            a = 0.0;
            b = 0.0;
            if j > 1 % need bottom
                a = a + (g(j)-g(j-1)+1/2*(sigma(j)+sigma(j-1))*-phi(j-1));
                b = b + -1/2*(sigma(j)+sigma(j-1));
            end
            if j < Ny % need top
                a = a + (-g(j+1)+g(j)-1/2*(sigma(j+1)+sigma(j))*phi(j+1));
                b = b + -1/2*(sigma(j+1)+sigma(j));
            end
            dphi(j) = phi(j) - a / b;
            phi(j) = a / b;
        end
        
        phi = phi - min(phi);
    
        if norm(dphi) < tol
            fprintf('iter = %i, error = %e\n', counter, norm(dphi))
            break;
        end
        counter = counter + 1;
    end
    
    % Calculate the fluxes
    for k = 1:nions
        Jy(2:Ny,k) = -D(k)*((rho(2:Ny,k)-rho(1:Ny-1,k))/dy + 0.5*Z(k)*(rho(2:Ny,k) + ...
            rho(1:Ny-1,k)).*(phi(2:Ny)-phi(1:Ny-1))/dy);
    end
    
    % Update internal rho values and calculate currents
    Iy = zeros(Ny+1,1);
    for k = 1:nions
        Iy = Iy + Z(k) * Jy(:,k);
        rho_new(:,k) = rho(:,k) - dt * (Jy(2:Ny+1,k) - Jy(1:Ny,k))/dy;
    end
    
    % Top Wall
    rho1 = rho(Ny,3);
    rho2 = rho(Ny,4);
    del1 = rho_new(Ny,3) - rho(Ny,3);
    del2 = rho_new(Ny,4) - rho(Ny,4);
    B = rho1 + rho2 + del1 + del2;
    C = rho1 * del2 + rho2 * del1 + del1 * del2;
    flux_dt = 1/2*(sqrt(B^2 - 4*C) - B);
    rho_new(Ny,3) = rho_new(Ny,3) + flux_dt;
    rho_new(Ny,4) = rho_new(Ny,4) + flux_dt;
    
    % Bottom Wall
    rho1 = rho(1,3);
    rho2 = rho(1,4);
    del1 = rho_new(1,3) - rho(1,3);
    del2 = rho_new(1,4) - rho(1,4);
    B = rho1 + rho2 + del1 + del2;
    C = rho1 * del2 + rho2 * del1 + del1 * del2;
    flux_dt = 1/2*(sqrt(B^2 - 4*C) - B);
    rho_new(1,3) = rho_new(1,3) + flux_dt;
    rho_new(1,4) = rho_new(1,4) + flux_dt;
        

    
    
%     epsilon = 80.2*8.8542e-12; % permittivity of water in F/m
%     e = 1.602e-19; % Fundamental charge in C
%     kB = 1.3806e-23; % Boltzmann constant in J/K
%     T = 293; % Temperature in K
%     
%     lapPhi = phi;
%     for i = 2:Nx-1
%         for j = 2:Ny-1
%             lapPhi(i,j) = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/dx^2+(phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/dy^2;
%         end
%     end
    

    
    if mod(counter2,10) == 0
        subplot(3,2,1)
        plot(y,rho(:,1))
%         colorbar
%         contour(xx,yy,rho(:,:,1).', [0.3 0.4 0.5 0.6 0.7 0.8 0.9])
        subplot(3,2,2)
        plot(y,rho(:,2))
%         contour(xx,yy,rho(:,:,2).', [0.2 0.4 0.6 0.8 1.0])
        subplot(3,2,3)
        plot(y,rho(:,3))
%         contour(xx,yy,rho(:,:,3).', [0.2 0.4 0.6 0.8 1.0])
        subplot(3,2,4)
        plot(y,rho(:,4))
%         contour(xx,yy,rho(:,:,4).', [0.2 0.4 0.6 0.8 1.0])
        subplot(3,2,5)
        plot(y,phi)
        subplot(3,2,6)
%         plot(y_Jy,Jy(:,1))
        plot(y,sumrhoz)
%         surf(x,y,-epsilon*kB*T/e^2/6.02e23*lapPhi.')
        pause(0.0000001)
    end
    counter2 = counter2 + 1;
    
    fprintf('time = %f, drho1/dt = %e\n', t, norm(rho_new(:,1) - rho(:,1)) / norm(rho(:,1)) / dt)
    rho = rho_new;
end
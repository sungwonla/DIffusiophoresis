clc;
clear;

L = 100e-6;
W = 100e-6;
% L = 1;
% W = 1;
nions = 4;
Nx = 21;
Ny = 21;

x = linspace(0,L,Nx);
y = linspace(0,W,Ny);
[xx,yy] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
x_Jx = linspace(0-dx/2,L+dx/2,Nx+1);
y_Jy = linspace(0-dy/2,W+dy/2,Ny+1);

% D = [1.330e-9 2.030e-9 9.310e-9 1.180e-9]; % Na, Cl, H, HCO3
D = [9.31e-9 2.03e-9 1.96e-9 2.01e-9]; % H+, Cl-, K+, Br-
% D = [9e-9 1e-9 9e-9 1e-9];
% D = D / max(D);
Z = [1 -1 1 -1];

dt = min(dx,dy)^2 / max(D) / 5;
tfinal = 10;
tol = 1e-9;

rho_H_high = 1;
rho_H_low = 0.01;

rho = zeros(Nx,Ny,nions);
phi = zeros(Nx,Ny);
temp = zeros(Nx,Ny);
Jx = zeros(Nx+1,Ny,nions);
Jy = zeros(Nx,Ny+1,nions);
Ix = zeros(Nx+1,Ny);
Iy = zeros(Nx,Ny+1);

w = 5e-6;
for i = 1:Nx
    for j = 1:Ny
        rho(i,j,1) = rho_H_low + (rho_H_high - rho_H_low) / L * x(i);
%         rho(i,j,1) = (rho_H_low+rho_H_high)/2 + (rho_H_high-rho_H_low)/2*erf((x(i)-L/2)/w);
        rho(i,j,3) = rho_H_low + (rho_H_high - rho_H_low) / L * y(j);
%         rho(i,j,3) = (rho_H_low+rho_H_high)/2 + (rho_H_high-rho_H_low)/2*erf((y(j)-L/2)/w);
    end
end
rho(:,:,2) = rho(:,:,1);
rho(:,:,4) = rho(:,:,3);

dphi = phi;
rho_new = rho;

counter2 = 1;
for t = 0:dt:tfinal

    g = zeros(Nx,Ny);
    sigma = zeros(Nx,Ny);
    sumrhoz = zeros(Nx,Ny);
    for k = 1:nions
        g = g + Z(k) * D(k) * rho(:,:,k);
        sigma = sigma + Z(k)^2 * D(k) * rho(:,:,k);
        sumrhoz = sumrhoz + Z(k) * rho(:,:,k);
    end
    
    % Find the electric potential by relaxation
    counter = 1;
    while true
        for i = 1:Nx
            for j = 1:Ny
                a = 0.0;
                b = 0.0;
                if i > 1 % need left
                    a = a + dy/dx*(g(i,j)-g(i-1,j)+1/2*(sigma(i,j)+sigma(i-1,j))*-phi(i-1,j));
                    b = b + -1/2*dy/dx*(sigma(i,j)+sigma(i-1,j));
                end
                if i < Nx % need right
                    a = a + dy/dx*(-g(i+1,j)+g(i,j)-1/2*(sigma(i+1,j)+sigma(i,j))*phi(i+1,j));
                    b = b + -1/2*dy/dx*(sigma(i+1,j)+sigma(i,j));
                end
                if j > 1 % need bottom
                    a = a + dx/dy*(g(i,j)-g(i,j-1)+1/2*(sigma(i,j)+sigma(i,j-1))*-phi(i,j-1));
                    b = b + -1/2*dx/dy*(sigma(i,j)+sigma(i,j-1));
                end
                if j < Ny % need top
                    a = a + dx/dy*(-g(i,j+1)+g(i,j)-1/2*(sigma(i,j+1)+sigma(i,j))*phi(i,j+1));
                    b = b + -1/2*dx/dy*(sigma(i,j+1)+sigma(i,j));
                end
                dphi(i,j) = phi(i,j) - a / b;
                phi(i,j) = a / b;
            end
        end
        
        phi = phi - min(min(phi));
        
        if norm(dphi) < tol
            fprintf('iter = %i, error = %e\n', counter, norm(dphi))
            break;
        end
        counter = counter + 1;
    end
    
    
    % Calculate the fluxes
%     phi = 0*phi;
    for k = 1:nions
        Jx(2:Nx,:,k) = -D(k)*((rho(2:Nx,:,k)-rho(1:Nx-1,:,k))/dx + 0.5*Z(k)*(rho(2:Nx,:,k) + ...
            rho(1:Nx-1,:,k)).*(phi(2:Nx,:)-phi(1:Nx-1,:))/dx);
        Jy(:,2:Ny,k) = -D(k)*((rho(:,2:Ny,k)-rho(:,1:Ny-1,k))/dy + 0.5*Z(k)*(rho(:,2:Ny,k) + ...
            rho(:,1:Ny-1,k)).*(phi(:,2:Ny)-phi(:,1:Ny-1))/dy);
    end
    
    % Update internal rho values
    drhodt = zeros(Nx,Ny,nions);
    for k = 1:nions
        rho_new(:,:,k) = rho(:,:,k) - dt * ((Jx(2:Nx+1,:,k) - Jx(1:Nx,:,k))/dx + ...
            (Jy(:,2:Ny+1,k) - Jy(:,1:Ny,k))/dy);
    end
    
    % Top Wall
    for i = 1:Nx
        rho1 = rho(i,Ny,3);
        rho2 = rho(i,Ny,4);
        del1 = rho_new(i,Ny,3) - rho(i,Ny,3);
        del2 = rho_new(i,Ny,4) - rho(i,Ny,4);
        B = rho1 + rho2 + del1 + del2;
        C = rho1 * del2 + rho2 * del1 + del1 * del2;
        flux_dt = 1/2*(sqrt(B^2 - 4*C) - B);
        rho_new(i,Ny,3) = rho_new(i,Ny,3) + flux_dt;
        rho_new(i,Ny,4) = rho_new(i,Ny,4) + flux_dt;
    end
    % Bottom Wall
    for i = 1:Nx
        rho1 = rho(i,1,3);
        rho2 = rho(i,1,4);
        del1 = rho_new(i,1,3) - rho(i,1,3);
        del2 = rho_new(i,1,4) - rho(i,1,4);
        B = rho1 + rho2 + del1 + del2;
        C = rho1 * del2 + rho2 * del1 + del1 * del2;
        flux_dt = 1/2*(sqrt(B^2 - 4*C) - B);
        rho_new(i,1,3) = rho_new(i,1,3) + flux_dt;
        rho_new(i,1,4) = rho_new(i,1,4) + flux_dt;
    end
    % Left Wall
    for j = 1:Ny
        rho1 = rho(1,j,1);
        rho2 = rho(1,j,2);
        del1 = rho_new(1,j,1) - rho(1,j,1);
        del2 = rho_new(1,j,2) - rho(1,j,2);
        B = rho1 + rho2 + del1 + del2;
        C = rho1 * del2 + rho2 * del1 + del1 * del2;
        flux_dt = 1/2*(sqrt(B^2 - 4*C) - B);
        rho_new(1,j,1) = rho_new(1,j,1) + flux_dt;
        rho_new(1,j,2) = rho_new(1,j,2) + flux_dt;
    end
    % Right Wall
    for j = 1:Ny
        rho1 = rho(Nx,j,1);
        rho2 = rho(Nx,j,2);
        del1 = rho_new(Nx,j,1) - rho(Nx,j,1);
        del2 = rho_new(Nx,j,2) - rho(Nx,j,2);
        B = rho1 + rho2 + del1 + del2;
        C = rho1 * del2 + rho2 * del1 + del1 * del2;
        flux_dt = 1/2*(sqrt(B^2 - 4*C) - B);
        rho_new(Nx,j,1) = rho_new(Nx,j,1) + flux_dt;
        rho_new(Nx,j,2) = rho_new(Nx,j,2) + flux_dt;
    end

    fprintf('time = %f, drho3/dt = %e\n', t, norm(rho_new(:,:,3) - rho(:,:,3)) / dt)
    rho = rho_new;
   
    Ix(2:Nx,:) = -(g(2:Nx,:)-g(1:Nx-1,:))/dx-1/2*(sigma(2:Nx,:)+sigma(1:Nx-1,:)).*(phi(2:Nx,:)-phi(1:Nx-1,:))/dx;
    Iy(:,2:Ny) = -(g(:,2:Nx)-g(:,1:Ny-1))/dy-1/2*(sigma(:,2:Nx)+sigma(:,1:Nx-1)).*(phi(:,2:Nx)-phi(:,1:Nx-1))/dy;
    Ix_alt = zeros(Nx+1,Ny);
    Iy_alt = zeros(Nx,Ny+1);
    for k = 1:nions
        Ix_alt = Ix_alt + Z(k)*Jx(:,:,k);
        Iy_alt = Iy_alt + Z(k)*Jy(:,:,k);
    end
    div_I = zeros(Nx,Ny);
    div_I_alt = zeros(Nx,Ny);
    for i = 1:Nx
        for j = 1:Ny
            div_I(i,j) = Ix(i+1,j)*dy-Ix(i,j)*dy+Iy(i,j+1)*dx-Iy(i,j)*dx;
            div_I_alt(i,j) = Ix_alt(i+1,j)*dy-Ix_alt(i,j)*dy+Iy_alt(i,j+1)*dx-Iy_alt(i,j)*dx;
        end
    end
    div_I = div_I / dx / dy;
    div_I_alt = div_I_alt / dx / dy;
    
    
    int_rho = [0 0 0 0];
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:nions
                int_rho(k) = int_rho(k) + rho(i,j,k) * dx * dy;
            end
        end
    end
    int_rho = int_rho / L / W;
    int_rho_mat(counter2,:) = int_rho;

    epsilon = 80.2*8.8542e-12; % permittivity of water in F/m
    e = 1.602e-19; % Fundamental charge in C
    kB = 1.3806e-23; % Boltzmann constant in J/K
    T = 293; % Temperature in K
    NA = 6.02e23; % Avogadro's number
    laplace_phi = zeros(Nx,Ny);
    dphidx = zeros(Nx,Ny);
    dphidy = zeros(Nx,Ny);

    for i = 2:Nx-1
        for j = 2:Ny-1
            laplace_phi(i,j) = (phi(i+1,j)-2*phi(i,j)+phi(i-1,j))/dx^2 + ...
                (phi(i,j+1)-2*phi(i,j)+phi(i,j-1))/dy^2;
            dphidx(i,j) = (phi(i+1,j)-phi(i-1,j))/2/dx;
            dphidy(i,j) = (phi(i,j+1)-phi(i,j-1))/2/dy;
        end
    end
    

    fprintf('int_rho = [%f %f %f %f]\n', int_rho)
    
%     pause(3)
    if mod(counter2,100) == 0
        subplot(3,2,1)
%         surf(x,y,rho(:,:,1).')
%         surf(x_Jx,y,Jx(:,:,4).')
        contour(xx,yy,rho(:,:,1).', [0.3 0.4 0.5 0.6 0.7 0.8 0.9])
        subplot(3,2,2)
%         surf(x,y,rho(:,:,2).')
        contour(xx,yy,rho(:,:,2).', [0.2 0.4 0.6 0.8 1.0])
        subplot(3,2,3)
%         surf(x,y_Jy,Jy(:,:,4).')
%         surf(x,y,rho(:,:,3).')
        contour(xx,yy,rho(:,:,3).', [0.2 0.4 0.6 0.8 1.0])
        subplot(3,2,4)
%         surf(x,y,rho(:,:,4).')
        contour(xx,yy,rho(:,:,4).', [0.2 0.4 0.6 0.8 1.0])
%         surf(x,y,-kB*T*epsilon/e^2/1000/6.02e23*laplace_phi.')
        subplot(3,2,5)
%         surf(x_Jx,y,Ix.')
        surf(x,y,phi.')
        subplot(3,2,6)
%         surf(x_Jx,y,Ix_alt.')
%         surf(x,y,drhodt(:,:,1).'-drhodt(:,:,2).'+drhodt(:,:,3).'-drhodt(:,:,4).')
%         surf(x,y,local_div.')
%         surf(x,y,div_I.')
%         surf(x,y,Z(1)*div_J1.'+Z(2)*div_J2.'+Z(3)*div_J3.'+Z(4)*div_J4.')
        surf(x,y,sumrhoz.')
%         surf(x,y,lapPhi.')
%         surf(x,y,rho(:,:,4).')
%         surf(x,y,-epsilon*kB*T/e^2/6.02e23*lapPhi.')
%         contour(xx,yy,phi.'*(T*kB/e)*1000, [10 20 30 40 50 60])
        axis square
        pause(3e-6)
    end
%     pause(3)
    counter2 = counter2 + 1;
%     return
end
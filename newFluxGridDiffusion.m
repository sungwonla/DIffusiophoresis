clc;
clf;
clear all;

%experimental parameters
d = [1.330e-9 2.030e-9 9.310e-9 1.180e-9]; %diffusivity constant: Na, Cl, H, HCO3
maxd = max(d); %finding maximum d to use in stability criterion
z = [1 -1 1 -1]; %charge on ion`
wallLengthX = 100e-6; %length of wall, x direction
wallLengthY = 100e-6; %length of wall, y direction
inletBulkVelocity = 0; %inlet velocity of bulk flow

%numerical approximation parameters
nodesX = 30; %number of nodes in x direction
nodesY = 30; %number of nodes in y direction
x = linspace(0,wallLengthX,nodesX); %x axis
y = linspace(0,wallLengthY,nodesY); %y axis
dx = x(2) - x(1); %distance step along x
dy = y(2) - y(1); %distance step along y
%dt = 0.5*((dx*dx)^2)/(((dx^2)+(dy^2))*maxd); %time step based on stability criterion
dt = dx^2/(5*maxd);
t_final = 1; %simulation runtime
tolerance = 1e-6;

%bulk velocity (one directional, along x)
bulkVelocityX = zeros(nodesY,nodesX); %bulk velocity values
bulkVelocityY = zeros(nodesY,nodesX); %flow is in one direction only
for j = 1:nodesY
    bulkVelocityX(j,:) = (-6*inletBulkVelocity/(wallLengthY^2))*(y(j)*(y(j) - wallLengthY));
end

rho = zeros(nodesY,nodesX,4); %rho values
phi = ones(nodesY,nodesX); %potential field values
fluxRhoX = zeros(nodesY,nodesX+1,4); %flux grid
fluxRhoY = zeros(nodesY+1,nodesX,4);
spaceCharge = zeros(nodesY,nodesX);
laplacianPotential = zeros(nodesY,nodesX); %laplacian of potential field

%initial conditions
minConcentration = 0.1;
maxConcentration = 1.0;
for i = 1:nodesX %concentration gradient along x
    rho(:,i,1) = minConcentration + (maxConcentration-minConcentration)*i/nodesX;
    rho(:,i,2) = minConcentration + (maxConcentration-minConcentration)*i/nodesX;
end

for j = 1:nodesY %concentration gradient along y
    rho(j,:,3) = minConcentration + (maxConcentration-minConcentration)*j/nodesY;
    rho(j,:,4) = minConcentration + (maxConcentration-minConcentration)*j/nodesY;
end


%updating mesh points
for t = 0:dt:t_final
    
    %node values for constraint equation
    sigma = ((z(1)^2)*d(1)*rho(:,:,1)+((z(2)^2)*d(2)*rho(:,:,2))+((z(3)^2)*d(3)*rho(:,:,3))+((z(4)^2)*d(4)*rho(:,:,4)));
    g = (z(1)*d(1)*rho(:,:,1))+(z(2)*d(2)*rho(:,:,2))+(z(3)*d(3)*rho(:,:,3))+(z(4)*d(4)*rho(:,:,4));
    A = (z(1)*rho(:,:,1))+(z(2)*rho(:,:,2))+(z(3)*rho(:,:,3))+(z(4)*rho(:,:,4));
    convergence = 1;
    
    %solving for potential with constraint equation
    while abs(convergence) > tolerance
        
        oldPhi = phi;
        
        %boundary points
        for i = 1:nodesX
            phi(nodesY,i) = phi(nodesY-1,i)+((g(nodesY,i)-g(nodesY-1,i))/sigma(nodesY,i)); %top wall
            phi(1,i) = phi(2,i)+((g(2,i)-g(1,i))/sigma(1,i)); %bottom wall
        end
        
        for j = 1:nodesY
            phi(j,nodesX) = phi(j,nodesX-1)-((g(j,nodesX)-g(j,nodesX-1))/sigma(j,nodesX))+(dx/sigma(j,nodesX))*(A(j,nodesX)*bulkVelocityX(j,nodesX)); %right wall
            phi(j,1) = phi(j,2)+((g(j,2)-g(j,1))/sigma(j,1))-((dx/sigma(j,1))*(A(j,1)*bulkVelocityX(j,1))); %left wall
        end
        
        %non-boundary points
        for i = 2:nodesX-1
            for j = 2:nodesY-1
                dsigmadx = differential(sigma(j,i-1),sigma(j,i+1),dx);
                dsigmady = differential(sigma(j-1,i),sigma(j+1,i),dy);
                dAdx = differential(A(j,i-1),A(j,i+1),dx);
                dAdy = differential(A(j-1,i),A(j+1,i),dy);
                laplaciang = laplacian(g(j,i),g(j,i-1),g(j,i+1),g(j-1,i),g(j+1,i),dx,dy);
                
                phi(j,i) = ((dsigmadx*phi(j,i+1)/dx)+(dsigmady*phi(j+1,i)/dy)+ ... %from gradient of phi
                    sigma(j,i)*(((phi(j,i+1)+phi(j,i-1))/(dx^2))+((phi(j+1,i)+phi(j-1,i))/(dy^2))) - ... %from laplacian of phi
                    dot([dAdx,dAdy],[bulkVelocityX(j,i),bulkVelocityY(j,i)]) + ... %gradA dot bulk velocity
                    laplaciang)/ ... %laplacian of g
                    ((dsigmadx/dx)+(dsigmady/dy)+(2*sigma(j,i)/(dx^2))+(2*sigma(j,i)/(dy^2))); %denominator term
            end
        end
        
        phi = phi - min(phi,[],'all');
        convergence = max(phi-oldPhi,[],'all'); %greatest change in potential value
        
    end
    
    %no-flux boundary conditions for flux grid
    fluxRhoX(:,1,:) = 0;
    fluxRhoX(:,nodesX+1,:) = 0;
    fluxRhoY(1,:,:) = 0;
    fluxRhoY(nodesY+1,:,:) = 0;
    
    %fluxes in x direction, inner
    for i = 2:nodesX
        for j = 1:nodesY
            for k = 1:4
                drhodx = (rho(j,i,k)-rho(j,i-1,k))/dx;
                dphidx = (phi(j,i)-phi(j,i-1))/dx;
                fluxRhoX(j,i,k) = -d(k)*drhodx - d(k)*z(k)*dphidx*0.5*(rho(j,i,k)+rho(j,i-1,k)) + 0.5*(rho(j,i,k)+rho(j,i-1,k))*0.5*(bulkVelocityX(j,i)+bulkVelocityX(j,i-1));
            end
        end
    end
    
    %fluxes in y direction, inner
    for i = 1:nodesX
        for j = 2:nodesY
            for k = 1:4
                drhody = (rho(j,i,k)-rho(j-1,i,k))/dy;
                dphidy = (phi(j,i)-phi(j-1,i))/dy;
                fluxRhoY(j,i,k) = -d(k)*drhody - d(k)*z(k)*dphidy*0.5*(rho(j,i,k)+rho(j-1,i,k));
            end
        end
    end
    
    %solving for rho using divergence theorem
    for i = 1:nodesX
        for j = 1:nodesY
            for k = 1:4
                rho(j,i,k) = rho(j,i,k) - (((fluxRhoX(j,i+1,k)-fluxRhoX(j,i,k))/dx) + ((fluxRhoY(j+1,i,k)-fluxRhoY(j,i,k))/dy))*dt;
            end
        end
    end
    
    %space charge
    for i = 1:nodesX
        for j = 1:nodesY
            spaceCharge(j,i) = z(1)*rho(j,i,1) + z(2)*rho(j,i,2) + z(3)*rho(j,i,3) + z(4)*rho(j,i,4);
        end
    end
    
    %laplacian of potential, corner
    laplacianPotential(1,1) = ((phi(1,3)-2*phi(1,2)+phi(1,1))/(dx^2)) + ((phi(3,1)-2*phi(2,1)+phi(1,1))/(dy^2));
    laplacianPotential(1,nodesX) = ((phi(1,nodesX-2)-2*phi(1,nodesX-1)+phi(1,nodesX))/(dx^2)) + ((phi(3,nodesX)-2*phi(2,nodesX)+phi(1,nodesX))/(dy^2));
    laplacianPotential(nodesY,1) = ((phi(nodesY,3)-2*phi(nodesY,2)+phi(nodesY,1))/(dx^2)) + ((phi(nodesY-2,1)-2*phi(nodesY-1,1)+phi(nodesY,1))/(dy^2));
    laplacianPotential(nodesY,nodesX) = ((phi(nodesY,nodesX-2)-2*phi(nodesY,nodesX-1)+phi(nodesY,nodesX))/(dx^2)) + ((phi(nodesY-2,nodesX)-2*phi(nodesY-1,nodesX)+phi(nodesY,nodesX))/(dy^2));
    %laplacian of potential, outer
    for i = 2:nodesX-1 %top and bottom
        laplacianPotential(1,i) = ((phi(1,i+1)-2*phi(1,i)+phi(1,i-1))/(dx^2)) + ((phi(3,i)-2*phi(2,i)+phi(1,i))/(dy^2));
        laplacianPotential(nodesY,i) = ((phi(nodesY,i+1)-2*phi(nodesY,i)+phi(nodesY,i-1))/(dx^2)) + ((phi(nodesY,i)-2*phi(nodesY-1,i)+phi(nodesY-2,i))/(dy^2));
    end
    for j = 2:nodesY-1 %left and right
        laplacianPotential(j,1) = ((phi(j+1,1)-2*phi(j,1)+phi(j-1,1))/(dx^2)) + ((phi(j,3)-2*phi(j,2)+phi(j,1))/(dy^2));
        laplacianPotential(j,nodesX) = ((phi(j+1,nodesX)-2*phi(j,nodesX)+phi(j-1,nodesX))/(dx^2)) + ((phi(j,nodesX)-2*phi(j,nodesX-1)+phi(j,nodesX-2))/(dy^2));
    end      
    %laplacian of potential, inner
    for i = 2:nodesX-1
        for j = 2:nodesY-1
            laplacianPotential(j,i) = laplacian(phi(j,i),phi(j,i-1),phi(j,i+1),phi(j-1,i),phi(j+1,i),dx,dy);
        end
    end
    
    
    %visualization
    mesh(x,y,phi);
    axis([0 wallLengthX 0 wallLengthY]);
    xlabel('Wall (Bottom)');
    ylabel('Wall (Left)');
    zlabel('rho');
    title(sprintf('Time = %f seconds',t));
    colorbar
    
    %fprintf('convergence = %f\n', convergence);
    fprintf('space charge = %f\n', min(spaceCharge,[],'all'))
    fprintf('laplacian potential = %f\n', min(laplacianPotential,[],'all'))
    
    pause(0.01);
    
end

function discreteDifferential = differential(fminus,fplus,d)
%gets numerical first derivative to scalar field
%fplus = f(j,i+1), fminus = f(j,i-1), d = dx or dy
    discreteDifferential = (fplus-fminus)/(2*d);
end

function discreteLaplacian = laplacian(f,fiminus,fiplus,fjminus,fjplus,dx,dy)
%gets laplacian of a scalar field
%fiminus = f(j,i-1), fiplus = f(j,i+1), fjminus = f(j-1,i), fjplus = f(j+1,i)
    discreteLaplacian = ((fiminus-(2*f)+fiplus)/(dx^2)) + ...
    ((fjminus-(2*f)+fjplus)/(dy^2));
end
clc;
clf;
clear all;

%experimental parameters
d1 = 1.330e-9; %diffusivity constant: Na, Cl, H, HCO3
d2 = 2.030e-9;
d3 = 9.310e-9;
d4 = 1.180e-9;
% d1 = 1.330; %diffusivity constant: Na, Cl, H, HCO3
% d2 = 2.030;
% d3 = 9.310;
% d4 = 1.180;
maxd = max([d1,d2,d3,d4]); %finding maximum d to use in stability criterion
z1 = 1; %charge on ion
z2 = -1;
z3 = 1;
z4 = -1;
wallLengthX = 1e-6; %length of wall, x direction
wallLengthY = 1e-6; %length of wall, y direction
inletBulkVelocity = 0; %inlet velocity of bulk flow

%numerical approximation parameters
nodesX = 80; %number of nodes in x direction
nodesY = 80; %number of nodes in y direction
x = linspace(0,wallLengthX,nodesX); %x axis
y = linspace(0,wallLengthY,nodesY); %y axis
dx = x(2) - x(1); %distance step along x
dy = y(2) - y(1); %distance step along y
dt = 0.5*((dx*dx)^2)/(((dx^2)+(dy^2))*maxd); %time step based on stability criterion
t_final = 0.3; %simulation runtime
tolerance = 0.0003;

%bulk velocity (one directional, along x)
bulkVelocityX = zeros(nodesY,nodesX); %bulk velocity values
bulkVelocityY = zeros(nodesY,nodesX); %flow is in one direction only
for j = 1:nodesY
    bulkVelocityX(j,:) = (-6*inletBulkVelocity/(wallLengthY^2))*(y(j)*(y(j) - wallLengthY));
end

rho1 = zeros(nodesY,nodesX); %rho values
drho1dt = zeros(nodesY,nodesX); %differential rho values
rho2 = zeros(nodesY,nodesX);
drho2dt = zeros(nodesY,nodesX);
rho3 = zeros(nodesY,nodesX);
drho3dt = zeros(nodesY,nodesX);
rho4 = zeros(nodesY,nodesX);
drho4dt = zeros(nodesY,nodesX);
phi = ones(nodesY,nodesX); %potential field values

%initial conditions
minConcentration = 0.1;
maxConcentration = 1.0;
for i = 1:nodesX %concentration gradient along x
    rho1(:,i) = minConcentration + (maxConcentration-minConcentration)*i/nodesX;
    rho2(:,i) = minConcentration + (maxConcentration-minConcentration)*i/nodesX;
end

for j = 1:nodesY %concentration gradient along y
    rho3(j,:) = minConcentration + (maxConcentration-minConcentration)*j/nodesY;
    rho4(j,:) = minConcentration + (maxConcentration-minConcentration)*j/nodesY;
end


%updating mesh points
for t = 0:dt:t_final
    
    %boundary conditions for rho values
%     for j = 1:nodesY
%         rho1(j,nodesX) = rho1(j,nodesX-1)/(1+(z1*(phi(j,nodesX)-phi(j,nodesX-1)))-(dx*bulkVelocityX(j,nodesX)/d1)); %right wall
%         rho1(j,1) = rho1(j,nodesX-1)/(1-(z1*(phi(j,1)-phi(j,2)))+(dx*bulkVelocityX(j,nodesX)/d1)); %left wall
%         
%         rho2(j,nodesX) = rho2(j,nodesX-1)/(1+(z2*(phi(j,nodesX)-phi(j,nodesX-1)))-(dx*bulkVelocityX(j,nodesX)/d2)); %right wall
%         rho2(j,1) = rho2(j,nodesX-1)/(1-(z2*(phi(j,1)-phi(j,2)))+(dx*bulkVelocityX(j,nodesX)/d2)); %left wall
%         
%         rho3(j,nodesX) = rho3(j,nodesX-1)/(1+(z3*(phi(j,nodesX)-phi(j,nodesX-1)))-(dx*bulkVelocityX(j,nodesX)/d3)); %right wall
%         rho3(j,1) = rho3(j,nodesX-1)/(1-(z3*(phi(j,1)-phi(j,2)))+(dx*bulkVelocityX(j,nodesX)/d3)); %left wall
%         
%         rho4(j,nodesX) = rho4(j,nodesX-1)/(1+(z4*(phi(j,nodesX)-phi(j,nodesX-1)))-(dx*bulkVelocityX(j,nodesX)/d4)); %right wall
%         rho4(j,1) = rho4(j,nodesX-1)/(1-(z4*(phi(j,1)-phi(j,2)))+(dx*bulkVelocityX(j,nodesX)/d4)); %left wall
%     end
%     
%     for i = 1:nodesX
%         rho1(nodesY,i) = rho1(nodesY-1,i)/(1+(z1*(phi(nodesY,i)-phi(nodesY-1,i)))); %top wall
%         rho1(1,i) = rho1(2,i)/(1-(z1*(phi(2,i)-phi(1,i)))); %bottom wall
%         
%         rho2(nodesY,i) = rho2(nodesY-1,i)/(1+(z2*(phi(nodesY,i)-phi(nodesY-1,i)))); %top wall
%         rho2(1,i) = rho2(2,i)/(1-(z2*(phi(2,i)-phi(1,i)))); %bottom wall
%         
%         rho3(nodesY,i) = rho3(nodesY-1,i)/(1+(z3*(phi(nodesY,i)-phi(nodesY-1,i)))); %top wall
%         rho3(1,i) = rho3(2,i)/(1-(z3*(phi(2,i)-phi(1,i)))); %bottom wall
%         
%         rho4(nodesY,i) = rho4(nodesY-1,i)/(1+(z4*(phi(nodesY,i)-phi(nodesY-1,i)))); %top wall
%         rho4(1,i) = rho4(2,i)/(1-(z4*(phi(2,i)-phi(1,i)))); %bottom wall
%     end
    
    
    %(fake) dirichlet boundary conditions
    rho1(nodesY,:) = 0.5; %top wall
    rho(1,:) = 0; %bottom wall
    rho1(:,1) = 0; %left wall
    rho1(:,nodesX) = 0; %right wall
   
    %neumann boundary conditions
    %rho(nodesY,:) = rho(nodesY-1,:); %top wall
    %rho1(1,:) = rho1(2,:); %bottom wall
    %rho(:,1) = rho(:,2); %left wall
    %rho(:,nodesX) = rho(:,nodesX-1) %right wall
    
    
    %node values for constraint equation
    sigma = ((z1^2)*d1*rho1)+((z2^2)*d2*rho2)+((z3^2)*d3*rho3)+((z4^2)*d4*rho4);
    g = (z1*d1*rho1)+(z2*d2*rho2)+(z3*d3*rho3)+(z4*d4*rho4);
    A = (z1*rho1)+(z2*rho2)+(z3*rho3)+(z4*rho4);
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
    
    
    %nernst-planck equations for non-boundary points
    for i = 2:nodesX-1
        for j = 2:nodesY-1
            
            dphidx = differential(phi(j,i-1),phi(j,i+1),dx);
            dphidy = differential(phi(j-1,i),phi(j+1,i),dy);
            laplacianPhi = laplacian(phi(j,i),phi(j,i-1),phi(j,i+1),phi(j-1,i),phi(j+1,i),dx,dy);
            
            %rho 1
            drho1dx = differential(rho1(j,i-1),rho1(j,i+1),dx);
            drho1dy = differential(rho1(j-1,i),rho1(j+1,i),dy);
            laplacianRho1 = laplacian(rho1(j,i),rho1(j,i-1),rho1(j,i+1),rho1(j-1,i),rho1(j+1,i),dx,dy);
            
            drho1dt = d1*laplacianRho1 - ... %diffusion term
                dot([bulkVelocityX(j,i),bulkVelocityY(j,i)],[drho1dx, drho1dy]) + ... %advection term
                d1*z1*dot([dphidx,dphidy],[drho1dx,drho1dy]) + ... %potential term 1
                d1*z1*rho1(j,i)*laplacianPhi; %potential term part 2
            rho1(j,i) = rho1(j,i) + drho1dt*dt;
            
            %rho 2
            drho2dx = differential(rho2(j,i-1),rho2(j,i+1),dx);
            drho2dy = differential(rho2(j-1,i),rho2(j+1,i),dy);
            laplacianRho2 = laplacian(rho2(j,i),rho2(j,i-1),rho2(j,i+1),rho2(j-1,i),rho2(j+1,i),dx,dy);
            
            drho2dt = d2*laplacianRho2 - ... %diffusion term
                dot([bulkVelocityX(j,i),bulkVelocityY(j,i)],[drho2dx, drho2dy]) + ... %advection term
                d2*z2*dot([dphidx,dphidy],[drho2dx,drho2dy]) + ... %potential term 1
                d2*z2*rho2(j,i)*laplacianPhi; %potential term part 2
            rho2(j,i) = rho2(j,i) + drho2dt*dt;
            
            %rho 3
            drho3dx = differential(rho3(j,i-1),rho3(j,i+1),dx);
            drho3dy = differential(rho3(j-1,i),rho3(j+1,i),dy);
            laplacianRho3 = laplacian(rho3(j,i),rho3(j,i-1),rho3(j,i+1),rho3(j-1,i),rho3(j+1,i),dx,dy);
            
            drho3dt = d3*laplacianRho3 - ... %diffusion term
                dot([bulkVelocityX(j,i),bulkVelocityY(j,i)],[drho3dx, drho3dy]) + ... %advection term
                d3*z3*dot([dphidx,dphidy],[drho3dx,drho3dy]) + ... %potential term 1
                d3*z3*rho3(j,i)*laplacianPhi; %potential term part 2
            rho3(j,i) = rho3(j,i) + drho3dt*dt;
            
            %rho 4
            drho4dx = differential(rho4(j,i-1),rho4(j,i+1),dx);
            drho4dy = differential(rho4(j-1,i),rho4(j+1,i),dy);
            laplacianRho4 = laplacian(rho4(j,i),rho4(j,i-1),rho4(j,i+1),rho4(j-1,i),rho4(j+1,i),dx,dy);
            
            drho4dt = d4*laplacianRho4 - ... %diffusion term
                dot([bulkVelocityX(j,i),bulkVelocityY(j,i)],[drho4dx, drho4dy]) + ... %advection term
                d4*z4*dot([dphidx,dphidy],[drho4dx,drho4dy]) + ... %potential term 1
                d4*z4*rho4(j,i)*laplacianPhi; %potential term part 2
            rho4(j,i) = rho4(j,i) + drho4dt*dt;
            
        end
    end
    
    %visualization
    mesh(x,y,rho2);
    axis([0 wallLengthX 0 wallLengthY]);
    xlabel('Wall (Bottom)');
    ylabel('Wall (Left)');
    zlabel('rho');
    title(sprintf('Time = %f seconds',t));
    colorbar
    
    fprintf('convergence = %f\n', t)
    
    %contour/2D color plots
    %     figure(2);
    %     subplot(3,1,1),contour(x,y,rho),colormap;
    %     subplot(3,1,2),pcolor(x,y,rho),shading interp,;
    %     subplot(3,1,3)';
    %     surf(rho');
    
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

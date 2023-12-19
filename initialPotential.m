clc;
clf;
clear all;

%experimental parameters
d1 = 1.330e-9; %diffusivity constant: Na, Cl, H, HCO3
d2 = 2.030e-9;
d3 = 9.310e-9;
d4 = 1.180e-9;
maxd = max([d1,d2,d3,d4]); %finding maximum d to use in stability criterion
z1 = 1; %charge on ion
z2 = -1;
z3 = 1;
z4 = -1;
wallLengthX = 1; %length of wall, x direction
wallLengthY = 1; %length of wall, y direction
inletBulkVelocity = 0; %inlet velocity of bulk flow

%numerical approximation parameters
nodesX = 30; %number of nodes in x direction
nodesY = 30; %number of nodes in y direction
x = linspace(0,wallLengthX,nodesX); %x axis
y = linspace(0,wallLengthY,nodesY); %y axis
dx = x(2) - x(1); %distance step along x
dy = y(2) - y(1); %distance step along y
dt = 0.5*((dx*dx)^2)/(((dx^2)+(dy^2))*maxd); %time step based on stability criterion
t_final = 0.3; %simulation runtime

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
potentialField = ones(nodesY,nodesX); %potential field values

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


%node values for constraint equation
sigma = ((z1^2)*d1*rho1)+((z2^2)*d2*rho2)+((z3^2)*d3*rho3)+((z4^2)*d4*rho4);
g = (z1*d1*rho1)+(z2*d2*rho2)+(z3*d3*rho3)+(z4*d4*rho4);
A = (z1*rho1)+(z2*rho2)+(z3*rho3)+(z4*rho4);
convergence = 1;

%solving for potential with constraint equation
while abs(convergence) > 0.0003
    oldPotentialField = potentialField;
    
    %boundary points
    for i = 1:nodesX
        potentialField(nodesY,i) = potentialField(nodesY-1,i)+((g(nodesY,i)-g(nodesY-1,i))/sigma(nodesY,i)); %top wall
        potentialField(1,i) = potentialField(2,i)+((g(2,i)-g(1,i))/sigma(1,i)); %bottom wall
    end
    
    for j = 1:nodesY
        potentialField(j,nodesX) = potentialField(j,nodesX-1)-((g(j,nodesX)-g(j,nodesX-1))/sigma(j,nodesX))+(dx/sigma(j,nodesX))*(A(j,nodesX)*bulkVelocityX(j,nodesX)); %right wall
        potentialField(j,1) = potentialField(j,2)+((g(j,2)-g(j,1))/sigma(j,1))-((dx/sigma(j,1))*(A(j,1)*bulkVelocityX(j,1))); %left wall
    end
    
    %non-boundary points
    for i = 2:nodesX-1
        for j = 2:nodesY-1
            dsigmadx = differential(sigma(j,i-1),sigma(j,i+1),dx);
            dsigmady = differential(sigma(j-1,i),sigma(j+1,i),dy);
            dAdx = differential(A(j,i-1),A(j,i+1),dx);
            dAdy = differential(A(j-1,i),A(j+1,i),dy);
            laplaciang = laplacian(g(j,i),g(j,i-1),g(j,i+1),g(j-1,i),g(j+1,i),dx,dy);
            
            potentialField(j,i) = ((dsigmadx*potentialField(j,i+1)/dx)+(dsigmady*potentialField(j+1,i)/dy)+ ... %from gradient of potential
                sigma(j,i)*(((potentialField(j,i+1)+potentialField(j,i-1))/(dx^2))+((potentialField(j+1,i)+potentialField(j-1,i))/(dy^2))) - ... %from laplacian of potential
                dot([dAdx,dAdy],[bulkVelocityX(j,i),bulkVelocityY(j,i)]) + ... %gradA dot bulk velocity
                laplaciang)/ ... %laplacian of g
                ((dsigmadx/dx)+(dsigmady/dy)+(2*sigma(j,i)/(dx^2))+(2*sigma(j,i)/(dy^2))); %denominator term
        end
    end
    
    potentialField = potentialField - min(potentialField,[],'all');
    
    convergence = max(potentialField-oldPotentialField,[],'all'); %greatest change in potential value
    fprintf('convergence = %f\n', convergence)
    
    %visualization
    mesh(x,y,potentialField);
    axis([0 wallLengthX 0 wallLengthY]);
    xlabel('Wall (Bottom)');
    ylabel('Wall (Left)');
    zlabel('potentialField');
    colorbar
    
    pause(0.01);
end

diffPhi = zeros(5,nodesX);
for i = 1:5
for j = 1:nodesX
    diffPhi(i,j) = (potentialField(i+1,j) - potentialField(i,j))/dx;
end
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

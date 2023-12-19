clc;
clear;

z = [1 -1 1 -1]; %charges
d = [1.330e-9 2.030e-9 9.310e-9 1.180e-9]; % Na, Cl, H, HCO3
e = 1.6e-19;
epsilon = 80;
kb = 1.38e-23;
T = 273;
endY = 100e-6; %end value of y
c1 = 0.001; %concentration at y=0
c2 = 1; %concentration at y=end

%guesses for fminsearch %dphi is 4.8e+05 or 2.1455e+05
grho1 = 6.15e-05; grho2 = 0.00924; grho3 = 0.0093; gdrho3 = 9750; gdphi = 2.1455e+05; %controllable variables
guesses = [grho1; grho2; grho3; gdrho3; gdphi];

%fminsearch finds minimized start boundary conditions
boundaryRhoStart = fminsearch(@(variableBoundaryRhos)runODE(z,e,epsilon,kb,T,c1,c2,d,endY,variableBoundaryRhos),guesses);
irho1 = boundaryRhoStart(1); irho2 = boundaryRhoStart(2); irho3 = boundaryRhoStart(3); 
idrho3 = boundaryRhoStart(4); idphi = boundaryRhoStart(5);

idrho1 = -irho1*z(1)*idphi; %dependent variables
idrho2 = -irho2*z(2)*idphi; 
irho4 = (c1*c1)/irho3;
idrho4 = ((-z(3)*d(3)/(z(4)*d(4)))*(idrho3+(irho3*z(3)*idphi))) - (irho4*z(4)*idphi);
iphi = 0;
boundaryRhoStart = [irho1,idrho1,irho2,idrho2,irho3,idrho3,irho4,idrho4,iphi,idphi];

%running ode45 with minimized start boundary conditions
[y,rho] = ode45(@(y,w) diff(y,w,z,e,epsilon,kb,T),[0,endY],boundaryRhoStart);
plot(y,rho(:,7));
xlabel('y'); ylabel('rho1');


function totalErr = runODE(z,e,epsilon,kb,T,c1,c2,d,endY,variableBoundaryRhos)
%returns total error at boundary y=end after running ODE45 (based on inputted B.Cs at y=0)

%variables to be run through fminsearch
vrho1 = variableBoundaryRhos(1);
vrho2 = variableBoundaryRhos(2);
vrho3 = variableBoundaryRhos(3); 
vdrho3 = variableBoundaryRhos(4);
vdphi = variableBoundaryRhos(5);

vdrho1 = -vrho1*z(1)*vdphi; %dependent variables %from no flux for Na
vdrho2 = -vrho2*z(2)*vdphi; %from no flux for Cl
vrho4 = (c1*c1)/vrho3; %from constant chemical potential for H, HCO3
vdrho4 = ((-z(3)*d(3)/(z(4)*d(4)))*(vdrho3+(vrho3*z(3)*vdphi))) - (vrho4*z(4)*vdphi); %from no current flux
vphi = 0;

boundaryRho = [vrho1; vdrho1; vrho2; vdrho2; vrho3; vdrho3; vrho4; vdrho4; vphi; vdphi];
[y,rho] = ode45(@(y,w) diff(y,w,z,e,epsilon,kb,T),[0,endY],boundaryRho);

%errors at boundary y=end
error1 = -d(1)*(rho(end,2)+(rho(end,1)*z(1)*rho(end,10))); %no flux for Na
error2 = -d(2)*(rho(end,4)+(rho(end,3)*z(2)*rho(end,10))); %no flux for Cl
error3 = rho(end,5)*rho(end,7) - (c2*c2); %constant chemical potential for H, HCO3
error4 = -z(3)*d(3)*(rho(end,6)+(rho(end,5)*z(3)*rho(end,10))) - (z(4)*d(4)*(rho(end,8)+(rho(end,7)*z(4)*rho(end,10)))); %no current flux
integralSum = 2*trapz([0 endY],[c1 c2]); %integral of Na + Cl of initial linear concetration profile
error5 = trapz(y,rho(:,1)) + trapz(y,rho(:,3)) - integralSum; %integral of Na + Cl
totalErr = abs(error1) + abs(error2) + abs(error3) + abs(error4) + abs(error5);
end


function dwdy = diff(y,w,z,e,epsilon,kb,T)
%diff eq, variables stored as w = [rho1, drho1dy, rho2, etc.]
rho1 = w(1); rho2 = w(3); rho3 = w(5); rho4 = w(7);
drho1dy = w(2); drho2dy = w(4); drho3dy = w(6); drho4dy = w(8);
phi = w(9); dphidy = w(10);

sumrho = z(1)*rho1 + z(2)*rho2 + z(3)*rho3 + z(4)*rho4;
dw1dy = drho1dy; dw3dy = drho2dy; dw5dy = drho3dy; dw7dy = drho4dy;
dw2dy = z(1)*(((rho1*(e^2)/(epsilon*kb*T))*sumrho) - drho1dy*dphidy);
dw4dy = z(2)*(((rho2*(e^2)/(epsilon*kb*T))*sumrho) - drho2dy*dphidy);
dw6dy = z(3)*(((rho3*(e^2)/(epsilon*kb*T))*sumrho) - drho3dy*dphidy);
dw8dy = z(4)*(((rho4*(e^2)/(epsilon*kb*T))*sumrho) - drho4dy*dphidy);
dw9dy = dphidy;
dw10dy = -((e^2)/(epsilon*kb*T))*sumrho;

dwdy = [dw1dy;dw2dy;dw3dy;dw4dy;dw5dy;dw6dy;dw7dy;dw8dy;dw9dy;dw10dy];
end
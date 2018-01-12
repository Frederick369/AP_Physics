%SPH 4U0
%Frederick Ngo, Lawrence Pang, Evan Li
%Mr. Van Bemmel
%Astrophysics Simulation, part 3

% Lay out all the constants
k = 1.38*10^-23 %JK^-1 Boltzmann constant
mu = 0.61
mh = 1.67*10^-27 %kg mass of hydrogen
%t = 1000000 % @ 0.8r, K temperature (initial value)
%rho = 800 %kgm^-3, mean density
a = 7.56*10^-16 % Jm^-3 K^-4 
G = 6.67*10^-11 % m^3kg^-1 s^-2
rs = 638000 % km, radius of sun-like object
cp = 2.5 %specific heat of pressure
cv = 1.5 %specific heat of volume
L = 3.86*10^26 % Luminosity
kappa =  0.001
c = 3*10^8 % ms^-1 speed of light
gamma = cp/cv

%Arbitrarily chosen dr
dr = 1000000
epsilon = 0.000442673749 % epsilon was reversed calculated
%p = k*rho*T/(mu*mh) % pressure
rho = ones(1,rs/dr); rho(1) = 86766.6819919466; % reverse calculated
press = ones(1,rs/dr); press(1) = (k/mu)/mh*rho(1); % pressure coefficients

%Note that the differentials are all /dr but not stated because
% / is not allowed in naming
% Initial values are chosen based on the initial constants outlined in the
% document
% Several initial values are reversed calculated. Although there will be
% errors in rounding initial values, a simulation is meant to model; it
% does not need to be exact

r = linspace(1000000,rs,1000000);
M = ones(1,rs/dr); M(1) = 1.09*10^24;
dM = ones(1,rs/dr); dM(1) = 4*pi*dr;
% Temperature values from core reverse calculated
% T actual value, chosen from the smaller of T calculated by radiative
% heating (Tr) and convective heating (Tc)
T = ones(1,rs/dr); T(1) = 2.64*10^6; 
Tr = ones(1,rs/dr); Tr(1) = 2.64*10^6;
Tc = ones(1,rs/dr); Tc(1) = 1.12*10^7;
dTr = ones(1,rs/dr); dTr(1) = -5.96*10^-4;
dTc = ones(1,rs/dr); dTc() = -2.15*10^-3;
dT = ones(1,rs/dr); dT(1) = min(dTr(1),dTc(1);
dP = ones(1,rs/dr); dP(1) = G*M(1)*rho(1)/dr^2;
P = ones(1,rs/dr); P(1) = G*L^2/r^4 + dP(1);

% Loop is started at 2 to prevent any out of bound errors when attempting
% to take previous values

for aa=2:r
    % Variable rho based on the approximation outlined in:
    % http://spacemath.gsfc.nasa.gov/Calculus/6Page102.pdf
    rho(aa)=e^(r(aa)^2*(-2.353338*10^-17)+10.38914695;
    dM(aa) = 4*pi*(r*aa)^2*rho(aa);
    M(aa) = M(aa-1) + dr*dM(aa);
    press(aa) = k/mu/mh*rho(aa);
    dP(aa) = G*M(aa)*rho(aa)/(aa*dr)^2;
    P(aa) = P(aa-1) + dr*dP(aa-1);
    dTr(aa) = (-3/4/a/c)*(kappa*rho(aa)/T(aa)^3)*(Lr(aa)/4/pi/r(aa)^2);
    dTc(aa) = (1-1/gamma)/press(aa)*dP;
    dT(aa) = min(dTr(aa),dTc(aa);
    Tr(aa) = Tr(aa-1) + dTr(aa-1)*dr;
    Tc(aa) = Tc(aa-1) + dTc(aa-1)*dr;       
end
% Superimposing other graphs onto this may be wanted, but this just plots
% the ideal temperature gradient
plot (dT,r)
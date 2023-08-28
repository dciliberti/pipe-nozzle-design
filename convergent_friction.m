% Preliminary design of an incompressible axisymmetric pipe/nozzle with friction

close all; clearvars; clc

% Input data
flowRate = 560;   % airflow in m3/hr
flowTemp = 40;      % initial flow temperature in Celsius
inletDiam = 0.17;   % inlet section diameter in m2
outletDiam = 0.13;  % outlet section diameter in m2
nozzleLen = 0.05;   % nozzle length in m
pressure = 101325;  % ambient pressure in Pa
gamma = 1.4;        % specific heat ratio
R = 287;            % gas constant in J/(kg K)
eps = 400;          % pipe/nozzle roughness height in micrometers
g = 9.81;           % gravitational accelearation in m/s2

% Anonymous functions
flowDensity = @(p,R,T) p / (R * (T + 273.15)); % flow density for a perfect gas
dynVisc = @(T) 1.458E-6 * (T+273.15)^1.5 / (T + 383.55); % Sutherland law for dynamic viscosity
reynoldsNum = @(rho,V,L,mu) rho*V*L/mu; % Reynolds number
areaRatio = @(M,gamma) (1/M) * ((2 + (gamma-1) * M^2) / (gamma+1)) ^ ((gamma+1)/(2*(gamma-1))); % area ratio for isentropic flow
% colebWhite = @(f,eps,D,Re) f^-0.5 + 2.0 * log10((eps / (3.71*D)) + (2.51 / (Re*sqrt(f))) ); % Colebrook and White formula for friction coefficient
haland = @(eps,D,Re) (-1.8 * log10((eps / (3.7*D))^1.11 + 6.9/Re))^-2; % Haland explicit formula of the Colebrook and White model

% Preliminary calculations
inletArea = pi*(inletDiam/2)^2;
outletArea = pi*(outletDiam/2)^2;
inletVel = flowRate/3600/inletArea;   % flow speed at inlet section (no friction effect) in m/s
outletVel = flowRate/3600/outletArea;   % flow speed at outlet section without friction in m/s
Mach = outletVel / sqrt(gamma*R*flowTemp); % Mach number at outlet section
if Mach > 0.3
    warning('The assumption of incompressible flow is not valid')
end
aStar = outletArea / areaRatio(Mach,gamma); % critical area in m2
if outletArea < aStar
    warning('Nozzle outlet area is less than the critical area')
end
flowDens = flowDensity(pressure,R,flowTemp); % flow density in kg/m3
flowVisc = dynVisc(flowTemp); % flow dynamic viscosity
staticPressure = pressure - 0.5*flowDens*inletVel^2; % static pressure at inlet
headTotal = (staticPressure / (flowDens*g)) + (inletVel^2 / (2*g)); % total head at inlet section

% Iterative design
n = 10;                 % divide the pipe/nozzle in n segments
diam = zeros(n+1,1);    % and n+1 sections
area = zeros(n+1,1);
long = nozzleLen/n;
for i = 0:n
    diam(i+1) = inletDiam + (inletDiam-outletDiam)*(cos(pi*(i/n))-1)/2; % cosinusoidal sections distribution
    area(i+1) = pi*(diam(i+1)/2)^2;
end

difference = 999; % initial value for convergence evaluation
tol = 1e-2; % tolerance for convergence
c = 0; % counter
deltaP = zeros(n+1,1);
while difference > tol
    c = c + 1;
    for i = 0:n
        vel = outletVel * area(i+1)/outletArea; % flow speed by continuity equation
        Re = reynoldsNum(flowDens,vel,diam(i+1),flowVisc);
        % f = fzero(@(f) colebWhite(f,eps*1E-6,diam,Re),0.02);
        f = haland(eps*1E-6,diam(i+1),Re);
        deltaP(i+1) = f * (long / diam(i+1)) * (flowDens * vel^2) / 2;    % local pressure loss in Pa
    end
    pressureLoss = sum(deltaP);
    headLoss = pressureLoss / (flowDens * g);
    flowSpeed = outletVel * sqrt(1 - headLoss / headTotal); % corrected flow speed with friction losses
    difference = abs(flowSpeed - outletVel); % evaluate difference between two consecutive steps
    outletVel = flowSpeed; % update old speed to new speed and iterate
end

s = stackedplot((0:n)/n,[diam/2,deltaP,cumsum(deltaP)],...
    "DisplayLabels",["Radius, m","Pressure loss, Pa","Cum. press. loss, Pa"]);
s.AxesProperties(1).YLimits = [0,inf];
s.AxesProperties(2).YLimits = [0,inf];
s.AxesProperties(3).YLimits = [0,inf];
xlabel('Non-dimensional pipe / nozze length')
grid on

frictionlessFlowSpeed = Mach*sqrt(gamma*R*flowTemp); % recap flow speed without friction from outlet Mach number estimated at first

disp('Comparison between flow speed with smooth (slip-wall) and rough pipe/nozzle in m/s:')
disp([frictionlessFlowSpeed,flowSpeed])
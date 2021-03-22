% Assignment-3
% Matthieu Bourbeau
% 100975211

%% PART 1
% Modified the Monte-Carlo simulator from assignment1 PART 2 by adding 
% a voltage across the x-axis of the semiconductor crystal. This voltage 
% results in an electric field forming within the
% semiconductor. 
% Set the applied voltage to 0.1 V, but increased to 0.4 V to clearly see 
% the effects of the electric field.

global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons

particles = 1000; 
timeStep = 1000;
time = zeros(1, timeStep);
tau = 0.2e-12;                      % mean time between collisions (s)
T = 300;                            % temperature (K)
vth = sqrt(2*C.kb*T/C.m_n);         % thermal velocity
dt = (100e-9)/vth/100;
temperature = zeros(1, timeStep);
Pscat = 1-exp(-dt/tau); 
region_x = linspace(0, 200, 400)*10^(-9); 
region_y = linspace(0, 100, 400)*10^(-9);  

voltagex = 0.1; 
voltagey = 0;
Ex = voltagex/(200e-9);
Ey = voltagey/(100e-9);

E = sqrt(Ex^2+Ey^2);
fprintf('The electric field seen by the electrons is equal to: %.3d V/m\n',E);

F = E*(C.q_0);
fprintf('The force on each electron is equal to: %.3d N\n',F);

a = F/C.m_n;
fprintf('The acceleration on each electron is equal to: %.3d m/s^2\n',a);

voltagex = 0.4; 
voltagey = 0;
Ex = voltagex/(200e-9);
Ey = voltagey/(100e-9);

Px = zeros(particles, timeStep);
Py = zeros(particles, timeStep); 

for n = 1: particles
    Px(n, 1) = region_x(randi(400));
    Py(n, 1) = region_y(randi(400));
end

Vx = zeros(particles, timeStep); 
Vy = zeros(particles, timeStep);
accelx = zeros(particles, timeStep); 
accely = zeros(particles, timeStep);

MaxwellBoltzmannVdist = makedist('Normal', 'mu', 0, 'sigma', sqrt(C.kb*T/C.m_n));

for k = 1: particles
    Vx(k, :) = random(MaxwellBoltzmannVdist);
    Vy(k, :) = random(MaxwellBoltzmannVdist);
    accelx(k, :) = Ex * (-C.q_0/C.m_n);
    accely(k, :) = Ey * (-C.q_0/C.m_n);
end

avgV = sqrt(sum(Vx(:, 1).^2)/particles + sum(Vy(:, 1).^2/particles));

for j = 1: particles
    for w = 2: timeStep
        Vx(j, w) = Vx(j, w-1) + accelx(j, w-1) * dt;
        Vy(j, w) = Vy(j, w-1) + accely(j, w-1) * dt;
        if isnan(Px(j, w-1))
            if left == 1
                Px(j, w) = 0 + Vx(j, w)*dt;
            end
            if right == 1
                Px(j, w) = 200e-9 + Vx(j, w)*dt;
            end
        else
            Px(j, w) = Px(j, w-1) + Vx(j, w)*dt;
        end
        
        if Px(j, w) > 200e-9
            left = 1; 
            right = 0;
            Px(j, w) = NaN;
        end
        if Px(j, w) < 0
            left = 0;
            right = 1;
            Px(j, w) = NaN;
        end
        
        Py(j, w) = Py(j, w-1) + Vy(j, w)*dt; 
        if Py(j, w) > 100e-9
            Py(j, w) = 100e-9; 
            Vy(j, w:end) = -Vy(j, w);
        end
        if Py(j, w) < 0
            Py(j, w) = 0;
            Vy(j, w:end) = -Vy(j, w);
        end
        if Pscat > rand()
            Vx(j, w:end) = random(MaxwellBoltzmannVdist);
            Vy(j, w:end) = random(MaxwellBoltzmannVdist);
        end
    end
end

for i = 1: timeStep
        temperature(i) = (sum(Vx(:, i).^2) + sum(Vy(:, i).^2))*C.m_n/C.kb/2/particles;
        if i > 1
            time(i) = time(i-1) + dt;
        end
end

n = 100e15;
for i = 1: timeStep
        Jd(i) = sqrt(sum(Vx(:, i).^2)/particles + sum(Vy(:, i).^2/particles))*(C.q_0);
        if i > 1
            time(i) = time(i-1) + dt;
        end
end

for g = 1: 10
    figure(1)
    plot(Px(g, :),Py(g, :))
    xlabel('x')
    ylabel('y')
    title('Plot of 2D Particle Trajectories')
    hold on
end

figure(2)
plot(time, Jd)
title('Drift Current vs. Time')
xlabel('Time (s)')
ylabel('Current (A)')

densityM = [Px(:, 1000), Py(:, 1000)];
figure(3)
hist3(densityM, [200 100])
title('Electron Density Map')
xlabel('x')
ylabel('y')

tempx = zeros(ceil(200),ceil(100));
tempy = zeros(ceil(200),ceil(100));
tempn = zeros(ceil(200),ceil(100)); 

for z = 1: particles
    x = floor(Px(z, 1000)/1e-9);
    y = floor(Py(z, 1000)/1e-9);
    if (x == 0 || isnan(x)) 
        x = 1; 
    end 
    if (y == 0 || isnan(y))
        y = 1; 
    end
    tempy(x, y) = tempy(x, y) + Vy(z, 1000)^2;
    tempx(x, y) = tempx(x, y) + Vx(z, 1000)^2;
    tempn(x, y) = tempn(x, y) + 1;
end

temp2 = (tempx + tempy).* C.m_n./C.kb./2./tempn;
temp2(isnan(temp2)) = 0; 
temp2 = temp2';

figure(4)
xtemp = linspace(1, 200, 200);
ytemp = linspace(1, 100, 100);
pcolor(xtemp, ytemp, temp2)
colormap(jet)
title('Temperature Density Map')
xlabel('x')
ylabel('y')
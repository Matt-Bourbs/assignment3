% Assignment-3
% Matthieu Bourbeau
% 100975211

%% PART 3
% Coupled PART 1 and PART 2.
% PART 2 was used to establish a potential of 0.8 V across the x-axis of the
% semiconductor and form the resulting electric field.
% PART 1 was modified to include the bottle neck barrier and determine
% the electrons' acceleration due to the electric field at each time step.

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
timeStep = 3500;
time = zeros(1, timeStep);
tau = 0.2e-12;                      % mean time between collisions (s)
T = 300;                            % temperature (K)
vth = sqrt(2*C.kb*T/C.m_n);         % thermal velocity
dt = (100e-9)/vth/500;
temperature = zeros(1, timeStep);
Pscat = 1-exp(-dt/tau); 
region_x = linspace(0, 200, 400)*10^(-9);
region_y = linspace(0, 100, 400)*10^(-9);

Px = zeros(particles, timeStep);
Py = zeros(particles, timeStep); 

for n = 1: particles
    Px(n, 1) = x(randi(400));
    Py(n, 1) = y(randi(400));
    while (Px(n, 1) >= 80e-9 && Px(n, 1) <= 120e-9 && Py(n, 1) >= 60e-9) || (Px(n, 1) >= 80e-9 && Px(n, 1) <= 120e-9 && Py(n, 1) <= 40e-9)
        Px(n, 1) = region_x(randi(400));
        Py(n, 1) = region_y(randi(400)); 
    end
end

Vx = zeros(particles, timeStep); 
Vy = zeros(particles, timeStep);
accelx = zeros(particles, timeStep); 
accely = zeros(particles, timeStep);

MaxwellBoltzmannVdist = makedist('Normal','mu',0,'sigma',sqrt(C.kb*T/C.m_n));

for k = 1: particles
    Vx(k, :) = random(MaxwellBoltzmannVdist);
    Vy(k, :) = random(MaxwellBoltzmannVdist);
    if round(Px(k, 1)*(10^9)) == 0 && round(Py(k, 1)*(10^9)) == 0 
        accelx(k, 1) = Ex(1, 1) * (-C.q_0/C.m_n); 
        accely(k, 1) = Ey(1, 1) * (-C.q_0/C.m_n);
    elseif round(Px(k, 1)*(10^9)) == 0 
        accelx(k, 1) = Ex(1, round(Py(k, 1)*(10^9))) * (-C.q_0/C.m_n); 
        accely(k, 1) = Ey(1, round(Py(k, 1)*(10^9))) * (-C.q_0/C.m_n);
    elseif round(Py(k, 1)*(10^9)) == 0
        accelx(k, 1) = Ex(round(Px(k, 1)*(10^9)), 1) * (-C.q_0/C.m_n); 
        accely(k, 1) = Ey(round(Px(k, 1)*(10^9)), 1) * (-C.q_0/C.m_n);
    else
        accelx(k, 1) = Ex(round(Px(k, 1)*(10^9)), round(Py(k, 1)*(10^9))) * (-C.q_0/C.m_n); 
        accely(k, 1) = Ey(round(Px(k, 1)*(10^9)), round(Py(k, 1)*(10^9))) * (-C.q_0/C.m_n);
    end
end

avgV = sqrt(sum(Vx(:, 1).^2)/particles + sum(Vy(:, 1).^2/particles));

for j = 1: particles
    for w = 2: timeStep
        Vx(j, w) = Vx(j, w-1) + accelx(j, w-1) * dt;
        Vy(j, w) = Vy(j, w-1) + accely(j, w-1) * dt;
        if isnan(Px(j, w-1))
            if left == 1
                if Vx(j, w) < 0
                    Vx(j, w:end) = -Vx(j, w);
                end
                Px(j, w) = 0 + Vx(j, w)*dt;
            end
            if right == 1
                if Vx(j, w) > 0
                    Vx(j, w:end) = -Vx(j, w);
                end
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
            
        if (Px(j, w) >= 80e-9 && Px(j, w) <= 120e-9 && Py(j, w) >= 60e-9) || (Px(j, w) >= 80e-9 && Px(j, w) <= 120e-9 && Py(j, w) <= 40e-9)
            if (Px(j, w-1) <= 80e-9 && Py(j, w-1) <= 40e-9) || (Px(j, w-1) <= 80e-9 && Py(j, w-1) >= 60e-9) || (Px(j, w-1) >= 120e-9 && Py(j, w-1) <= 40e-9) || (Px(j, w-1) >= 120e-9 && Py(j, w-1) >= 40e-9)
                Vx(j, w:end) = -Vx(j, w-1);
            end
            if Px(j, w-1) >= 80e-9 && Px(j, w-1) <= 120e-9 && Py(j, w-1) <= 60e-9 && Py(j, w-1) >= 40e-9
                Vy(j, w:end) = -Vy(j, w-1);
            end
        end
        
        if Pscat > rand()
            Vx(j, w:end) = random(MaxwellBoltzmannVdist);
            Vy(j, w:end) = random(MaxwellBoltzmannVdist);
        end
        
        if isnan(Px(j, w))
            if round(Px(j, w-1)*(10^9)) == 0 && round(Py(j, w-1)*(10^9)) == 0 
                accelx(j, w) = Ex(1, 1) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(1, 1) * (-C.q_0/C.m_n);
            elseif round(Px(j, w-1)*(10^9)) == 0 
                accelx(j, w) = Ex(1, round(Py(j, w-1)*(10^9))) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(1, round(Py(j, w-1)*(10^9))) * (-C.q_0/C.m_n);
            elseif round(Py(j, w-1)*(10^9)) == 0
                accelx(j, w) = Ex(round(Px(j, w-1)*(10^9)), 1) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(round(Px(j, w-1)*(10^9)), 1) * (-C.q_0/C.m_n);
            else 
                accelx(j, w) = Ex(round(Px(j, w-1)*(10^9)), round(Py(j, w-1)*(10^9))) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(round(Px(j, w-1)*(10^9)), round(Py(j, w-1)*(10^9))) * (-C.q_0/C.m_n);
            end
        else
            if round(Px(j, w)*(10^9)) == 0 && round(Py(j, w)*(10^9)) == 0 
                accelx(j, w) = Ex(1, 1) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(1, 1) * (-C.q_0/C.m_n);
            elseif round(Px(j, w)*(10^9)) == 0 
                accelx(j, w) = Ex(1, round(Py(j, w)*(10^9))) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(1, round(Py(j, w)*(10^9))) * (-C.q_0/C.m_n);
            elseif round(Py(j, w)*(10^9)) == 0
                accelx(j, w) = Ex(round(Px(j, w)*(10^9)), 1) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(round(Px(j, w)*(10^9)), 1) * (-C.q_0/C.m_n);
            else 
                accelx(j, w) = Ex(round(Px(j, w)*(10^9)), round(Py(j, w)*(10^9))) * (-C.q_0/C.m_n); 
                accely(j, w) = Ey(round(Px(j, w)*(10^9)), round(Py(j, w)*(10^9))) * (-C.q_0/C.m_n);
            end
        end
    end
end

for g = 1:10
    figure(1)
    xLim = [80e-9 80e-9 120e-9 120e-9]; 
    yLim1 = [0 40e-9 40e-9 0];
    yLim2 = [100e-9 60e-9 60e-9 100e-9];
    line(xLim, yLim1)
    hold on
    line(xLim, yLim2)
    plot(Px(g, :), Py(g, :))
    xlabel('x (nm)')
    ylabel('y (nm)')
    title('Plot of 2D Particle Trajectories')
    hold off
end

% Considering there is an electric field present in the semiconductor, the
% expected behaviour of the electrons is that they will tend to travel
% towards the right side of the crystal. This effect can be seen in the
% following plots showing a 3D density plot and a 2D density map:

densityM = [Px(:, 1000), Py(:, 1000)];
figure(2)
hist3(densityM, [200 100])
title('Electron Density Map')
xlabel('x')
ylabel('y')
% Assignment-3
% Matthieu Bourbeau
% 100975211

%% PART 2
% Modified assigment2 to calculate and plot the potential over the 
% semiconductor crystal with a bottle-neck barrier present. 
% Set the potential across the x-axis to 0.8 V to clearly demonstrate 
% the effects. 
% The electric field was calculated from the potential and plotted. 

L = 200; 
W = 100; 
dx = 1; 
dy = 1; 
nx = L/dx; 
ny = W/dy;

G = sparse(nx*ny, nx*ny);
F = zeros(1, nx*ny); 

sigma = ones(nx, ny); 
for i = 1:nx
    for j = 1:ny
        if j <= (40) || j >= (60)
            if i >= (80) && i <= (120)
                sigma(i, j) = 10^(-2); 
            end
        end
    end
end

for i = 1: nx
    for j = 1:ny       
        n = j + (i - 1) * ny;         
        if i == 1
            G(n, :) = 0;
            G(n, n) = 1; 
            F(n) = 0.8;
        elseif i == nx
            G(n, n) = 1; 
        elseif j == 1
            sigmaUPPER = (sigma(i, j) + sigma(i, j+1)) / 2.0;
            sigmaRIGHT = (sigma(i, j) + sigma(i+1, j)) / 2.0;
            sigmaLEFT = (sigma(i, j) + sigma(i-1, j)) / 2.0;
            
            G(n, n) = -sigmaUPPER - sigmaRIGHT - sigmaLEFT; 
            G(n, n + 1) = sigmaUPPER;
            G(n, n + ny) = sigmaRIGHT;
            G(n, n - ny) = sigmaLEFT;
        elseif j == ny
            sigmaRIGHT = (sigma(i, j) + sigma(i+1, j)) / 2.0;
            sigmaLEFT = (sigma(i, j) + sigma(i-1, j)) / 2.0;
            sigmaLOWER = (sigma(i, j) + sigma(i, j-1)) / 2.0; 
            
            G(n, n) = -sigmaUPPER - sigmaRIGHT - sigmaLEFT; 
            G(n, n - 1) = sigmaLOWER;
            G(n, n + ny)= sigmaRIGHT;
            G(n, n - ny)= sigmaLEFT;
        else
            sigmaUPPER = (sigma(i, j) + sigma(i, j+1)) / 2.0;
            sigmaRIGHT = (sigma(i, j) + sigma(i+1, j)) / 2.0;
            sigmaLEFT = (sigma(i, j) + sigma(i-1, j)) / 2.0;
            sigmaLOWER = (sigma(i, j) + sigma(i, j-1)) / 2.0; 
            
            G(n, n)= -sigmaUPPER - sigmaLOWER - sigmaRIGHT - sigmaLEFT;
            G(n, n + 1) = sigmaUPPER;
            G(n, n - 1) = sigmaLOWER;
            G(n, n + ny) = sigmaRIGHT;
            G(n, n - ny) = sigmaLEFT;
        end
    end
end

V2 = zeros(nx, ny); 
V = G\F'; 

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1)*ny;
        V2(i, j) = V(n); 
    end
end

for i = 1:nx
    for j = 1:ny
        if sigma(i, j) == power(10, -12)
            V2(i, j) = sigma(i, j); 
        end
    end
end

figure(1)
mesh(V2)
colormap(jet)
title('Potential with Bottle-Neck Region')
xlabel('y')
ylabel('x')

[Ey, Ex] = gradient(V2*1e9); 
Ex = -Ex; 
Ey = -Ey; 

figure(2)
quiver(Ex, Ey)
camroll(90)
xlim([0 100])
ylim([0 200])
title('Electric Field')
xlabel('y')
ylabel('x')
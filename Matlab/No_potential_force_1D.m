clear all; close all;

%% Parameters
a = 1;  % Length
M = 1; % Mass
N = 100; % Spatial grid points
dt = 1e-3; % Time step
NPt = 5000; % Number of time steps
V0 = 1; % Potential barrier height
hbar = 1;  % Planck's constant
m = 1;      % Mass of particle

%% Spatial and Momentum Grids
x = linspace(-a/2, a/2, N)'; % Spatial grid
dx = x(2) - x(1);        % Spatial step size
k = fftshift(2 * pi * (0:N-1)' / (N * dx)) - pi/dx; % Centered momentum grid

%% Time Range
T = dt * NPt;
t = linspace(0, T, NPt);

%% Force Parameters
F0 = 1; % Force amplitude
gamma = 0.05; % Damping factor for the Gaussian envelope of the force
t0 = 2;

% Define the time-dependent force function
F_t = @(t) F0 * exp(-gamma * t) .* (t >= t0);

%% Initial Wave Packet
width = a / 10; % Width parameter for the Gaussian spread
x0 = 0; % Center of the wave packet

Phi0 = exp(-((x - x0).^2) / (2 * width^2)); % Gaussian wave packet
Phi0 = exp(-(5*(x-0*a/128)).^2); 
%% Kinetic and Potential Operators for Split-Step Method
kinetic_energy = (hbar^2 * k.^2) / (2 * m);
GK_half = fftshift(exp(-1i * (dt / (2 * hbar)) * kinetic_energy)); % Half-step kinetic propagator
GK_full = fftshift(exp(-1i * (dt / hbar) * kinetic_energy)); % Full-step kinetic propagator

%% Static Potential (Barrier) and Force-dependent Potential
V = -V0 * ones(N, 1); % Constant potential background
V0 = 200;
V = zeros(length(x),1) - V0; % 1*((2*x).^2 - (0.6*a)^2);  %  
b = a/16;
V(x<-b) = 0;
V(x>+b) = 0;

GV = exp(-0.5i * (dt / hbar) * V); % Static potential operator

%% Initialize Wave Packet
Phi = Phi0;
Prob_density = zeros(N, NPt); % Store probability density evolution

%% Time Evolution with Split Operator Method
for nrn = 1:NPt
    % Time-dependent force term in real space
    force_term = F_t(t(nrn)) * x;
    GV_force = exp(-1i * (dt / (2 * M)) * force_term / hbar); % Time-dependent force

    % Step 1: Half-step potential in real space
    Phi = GV .* GV_force .* Phi;

    % Step 2: Full-step kinetic evolution in momentum space
    Phi_k = fft(Phi); 
    Phi_k = GK_full .* Phi_k;  % Apply full-step kinetic propagator
    Phi = ifft(Phi_k); 

    % Step 3: Half-step potential in real space
    Phi = GV .* GV_force .* Phi;

    % Normalize wavefunction
    norm_factor = sqrt(sum(abs(Phi).^2) * dx);
    Phi = Phi / norm_factor;
    % Store probability density at each time step
    Prob_density(:, nrn) = abs(Phi).^2;
    % Plot evolution of probability density at each time step
    if mod(nrn, 10) == 0 % Plot every 10th time step
        plot(x, abs(Phi).^2, 'LineWidth', 2);
        xlabel('Position (x)');
        ylabel('Probability Density |Φ(x,t)|^2');
        title(['Time Evolution of Wave Packet at Time ', num2str(t(nrn))]);
        axis([-a/2 a/2 0 max(abs(Phi0).^2)]);
        grid on;
        pause(0.1); % Pause to create animation effect
    end
end

%% Plot Probability Density Evolution
figure;
plot(x, abs(Phi0).^2, 'LineWidth', 2); % Initial wave packet
hold on;
plot(x, abs(Phi).^2, 'r--', 'LineWidth', 2); % Final wave packet
xlabel('Position (x)');
ylabel('Probability Density |Φ(x)|^2');
title('Initial and Final Wave Packet');
legend('Initial Wave Packet', 'Final Wave Packet');
grid on;

%% 3D Surface Plot of Probability Density
[X, T] = meshgrid(x, t);
figure;
surf(T', X', Prob_density, 'EdgeColor', 'none');
xlabel('Time');
ylabel('Position (x)');
zlabel('Probability Density |Ψ(x, t)|^2');
title('Probability Density Evolution Over Time');
colorbar;
view(3);
shading interp;

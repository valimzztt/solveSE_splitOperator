clear all; close all;

%% Parameters
a = 1;  % Length
M = 1; % Mass
NPt = 5000; % Number of time steps
V0 = 1; % Potential barrier height
hbar = 1;  % Planck's constant
m = 1;      % Mass of particle
% Constants (Metric Units)
N = 512;          % Number of points to use
L = 4e-9;         % Extent of simulation (in meters)
%% Spatial and Momentum Grids
x =  L * linspace(-0.5, 0.5 - 1.0 / N, N)';  % Spatial grid
dx = x(2) - x(1);        % Spatial step size
k = fftshift(2 * pi * (0:N-1)' / (N * dx)) - pi/dx; % Centered momentum grid

%% Time Range
dt = 1e-3;       % Timestep in seconds
T = dt * NPt;
t = linspace(0, T, NPt);



%% Force Parameters
F0 = 0; % Force amplitude
gamma = 0.05; % Damping factor for the Gaussian envelope of the force
t0 = 2;

% Define the time-dependent force function
F_t = @(t) F0 * exp(-gamma * t) .* (t >= t0);

%% Initial Wave Packet
% Initial Wavefunction (Gaussian)
SIGMA = 0.07;
x0 = -0.0;
Phi0 = exp(-((x / L - x0) / SIGMA).^2 / 2.0);
Phi0 = Phi0 / sqrt(sum(abs(Phi0.^2) * dx));  % Normalize wavefunction
%% Kinetic and Potential Operators for Split-Step Method
kinetic_energy = (hbar^2 * k.^2) / (2 * m);
GK_half = fftshift(exp(-1i * (dt / (2 * hbar)) * kinetic_energy)); % Half-step kinetic propagator
GK_full = fftshift(exp(-1i * (dt / hbar) * kinetic_energy)); % Full-step kinetic propagator

%% Static Potential (HO) and Force-dependent Potential
% Uncomment one of the following potentials as needed

% Finite Square Well Potential
% V = ones(N, 1);
% V(abs(X) <= L / 4.0) = 0.0;

% Simple Harmonic Oscillator Potential
V = 15  * (x / L).^2;

% Alternative Potentials (Uncomment as needed)
% V = 3.5 * 1e-18 * (X / L); % Linear Potential
% V = 15 * 1e-18 * ((X / L).^2 + exp(-0.5 * (X / L).^2 / 0.05^2) / 8.0); % Complex Potential

% Plot Initial Wavefunction and Potential
figure;
subplot(2,1,1);
plot(x * 1e9, abs(Phi0).^2, 'LineWidth', 2);
xlabel('Position (nm)');
ylabel('Probability Density |Ψ(x)|^2');
title('Initial Wavefunction');
grid on;


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

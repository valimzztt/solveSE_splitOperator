clear all; close all;

%% Parameters
a = 1;             % Length
M = 1;             % Mass
NPt = 5000;        % Number of time steps
V0 = 1;            % Potential barrier height
hbar = 1;          % Planck's constant
m = 1;             % Mass of particle
N = 512;           % Number of points
L = 4e-9;          % Extent of simulation (meters)

%% Spatial and Momentum Grids
x = L * linspace(-0.5, 0.5 - 1.0 / N, N)';  % Spatial grid
dx = x(2) - x(1);                           % Spatial step size
k = (2 * pi / (N * dx)) * [0:N/2-1 -N/2:-1]'; % Momentum grid

%% Time Range
dt = 1e-3;        % Timestep
T = dt * NPt;
t = linspace(0, T, NPt);

%% Initial Wave Packet (Gaussian)
SIGMA = 0.07;
x0 = -0.0;
Phi0 = exp(-((x / L - x0) / SIGMA).^2 / 2.0);
Phi0 = Phi0 / sqrt(sum(abs(Phi0).^2) * dx);  % Normalize

%% Potential (Harmonic Oscillator as Example)
V = 0.5 * (15 * 1e-18 * (x / L).^2); % Adjust scale if needed

%% Define Time-Dependent Force
F0 = 0;       % Force amplitude
gamma = 0.05; % Damping factor
t0 = 2;
F_t = @(t) F0 * exp(-gamma * t) .* (t >= t0);

%% Split-Step Propagators
kinetic_energy = (hbar^2 * k.^2) / (2 * m);
GK_half = exp(-1i * (dt / (2 * hbar)) * kinetic_energy);  % Half-step kinetic
GK_full = exp(-1i * (dt / hbar) * kinetic_energy);         % Full-step kinetic
GV = exp(-0.5i * (dt / hbar) * V);                         % Static potential

%% Initialize Wave Packet
Phi = Phi0;
Prob_density = zeros(N, NPt); % Store probability density evolution

%% Time Evolution
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
    assert(abs(norm_factor - 1) < 1e-6, 'Normalization deviates at step %d', nrn);

    % Store probability density
    Prob_density(:, nrn) = abs(Phi).^2;
     % Plot evolution of probability density at each time step
    if mod(nrn, 100) == 0 % Plot every 100th time step
        plot(x, abs(Phi).^2, 'LineWidth', 2);
        xlabel('Position (x)');
        ylabel('Probability Density |Φ(x,t)|^2');
        title(['Time Evolution of Wave Packet at Time ', num2str(t(nrn))]);
        axis([-a/2 a/2 0 max(abs(Phi0).^2)]);
        grid on;
        pause(0.1); % Pause to create animation effect
    end
end

%% Plot Initial and Final Wave Packets
figure;
plot(x, abs(Phi0).^2, 'LineWidth', 2); % Initial wave packet
hold on;
plot(x, abs(Phi).^2, 'r--', 'LineWidth', 2); % Final wave packet
xlabel('Position (x)');
ylabel('Probability Density |Φ(x)|^2');
title('Initial and Final Wave Packet');
legend('Initial Wave Packet', 'Final Wave Packet');
grid on;

%% 3D Surface Plot of Probability Density Evolution
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

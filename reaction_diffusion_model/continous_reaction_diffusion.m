clear all; close all; clc;

% --- Standard Paper Constants ---
C = 1.35;
beta = 4.0;
% ECG Weights: Linear combination from Equation (4)
alpha = [-0.024, 0.0216, -0.0012, 0.12];

% Initial Conditions (Small kick to start)
x0 = [0.0; 0.0; 0.1; 0.0];

% 1. Simulation Settings (Run longer to eliminate startup jump)
% [0 to 300] ensures we have stable chaos for VF
tspan = [0 300];

% 2. Numerical Precision Options (MUST have for H=2.164)
% Tightening RelTol to 1e-8 prevents the flatline collapse
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% --- Simulation & Plotting ---
H_values = [3.0, 2.729, 2.164];
subplot_order = [1, 2, 3];
panels = {'(a)', '(c)', '(e)'};
titles = {'Normal Rhythm', 'Quasiperiodicity', 'VF'};

figure('Name', 'Replication of Figure 4 (ECG panels)', 'Position', [100, 100, 1000, 800]);

HR_bpm = 90;
Gamma_t = 0.08804 * HR_bpm  - 0.06754;

for i = 1:3
    H = H_values(i);

    % Define the dynamic system in (x1, x2, x3, x4)
    % dx/dt as a column vector
    ode_fun = @(t, x) Gamma_t * [ ...
        x(1) - x(2) - C*x(1)*x(2) - x(1)*x(2)^2; ...                                  % dx1/dt
        H*x(1) - 3*x(2) + C*x(1)*x(2) + x(1)*x(2)^2 + beta*(x(4) - x(2)); ...         % dx2/dt
        x(3) - x(4) - C*x(3)*x(4) - x(3)*x(4)^2; ...                                  % dx3/dt
        H*x(3) - 3*x(4) + C*x(3)*x(4) + x(3)*x(4)^2 + 2*beta*(x(2) - x(4)) ...        % dx4/dt 
    ];

    % Solve the ODE system
    [t, x_out] = ode15s(ode_fun, tspan, x0, options);

    % 3. Apply Equation (4): Linear mixture of the four components
    ECG = alpha(1)*x_out(:,1) + alpha(2)*x_out(:,2) + alpha(3)*x_out(:,3) + alpha(4)*x_out(:,4);

    % 4. Create Subplots to match paper layout
    subplot(3, 1, subplot_order(i));

    % Zoom tightly into the steady state (Showing only ~5 seconds)
    idx_zoom = find(t > 195 & t < 200); % Matches paper x-axis exactly

    % Plot with a thin line as in the paper
    plot(t(idx_zoom), ECG(idx_zoom), 'k', 'LineWidth', 2.0);

    % Annotation and Labeling
    axis tight;
    title(['Figure 4', panels{i}, ' - ', titles{i}, ' (H = ', num2str(H), ')']);
    ylabel('\xi(t) [ECG]');
    grid on;
    if i == 3; xlabel('Time (t)'); end
end

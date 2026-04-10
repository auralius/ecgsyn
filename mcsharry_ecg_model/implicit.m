clear all; close all; clc;

% ============================================================
% ECGSYN - Fixed-step implicit Tustin with full analytic Jacobian
% ============================================================

% -----------------------------
% Simulation settings
% -----------------------------
dt = 1/256;                 % fixed step
t_end = 30;
t = 0:dt:t_end;
n_steps = length(t);

max_iter = 6;
tol = 1e-10;

% -----------------------------
% ECGSYN parameters
% -----------------------------
fhi = 0.25;                 % baseline wander frequency
A   = 0.005;                % baseline wander amplitude

% Morphology defaults from ECGSYN
ti_deg = [-60, -15, 0, 15, 90];
ai = [1.2, -5.0, 30.0, -7.5, 0.75];
bi = [0.25, 0.10, 0.10, 0.10, 0.40];

% Convert angles to radians
ti = ti_deg * pi / 180;

% Optional heart-rate scaling as in ECGSYN
hr_mean = 60.0;
hrfact  = sqrt(hr_mean / 60.0);
hrfact2 = sqrt(hrfact);

bi = bi * hrfact;
ti(1) = ti(1) * hrfact2;
ti(2) = ti(2) * hrfact;
ti(3) = ti(3) * 1.0;
ti(4) = ti(4) * hrfact;
ti(5) = ti(5) * 1.0;

% -----------------------------
% Omega profile
% Start with constant omega for clean solver comparison.
% Replace later with RR-driven omega(t) if desired.
% -----------------------------
omega0 = 2*pi;              % 1 beat per second
omega = omega0 * ones(1, n_steps);

% -----------------------------
% Initial condition
% -----------------------------
x = zeros(3, n_steps);      % [x; y; z]
x(:,1) = [1.0; 0.0; 0.04];

% ============================================================
% Main implicit Tustin loop
% ============================================================
for n = 1:n_steps-1
    xp = x(:, n);           % previous state
    xn = xp;                % initial guess for Newton

    tp = t(n);
    tn = t(n+1);

    omega_p = omega(n);
    omega_n = omega(n+1);

    f_p = ecgsyn_rhs(xp, tp, omega_p, ti, ai, bi, fhi, A);

    for iter = 1:max_iter
        f_n = ecgsyn_rhs(xn, tn, omega_n, ti, ai, bi, fhi, A);

        % Tustin residual
        R = xn - xp - (dt/2) * (f_n + f_p);

        % Full analytic Jacobian of residual wrt xn
        Jf = ecgsyn_rhs_jacobian(xn, tn, omega_n, ti, ai, bi, fhi, A);
        J  = eye(3) - (dt/2) * Jf;

        dx = J \ R;
        xn = xn - dx;

        if norm(dx, inf) < tol
            break;
        end
    end

    x(:, n+1) = xn;
end

z = x(3,:);

% ============================================================
% Plot
% ============================================================
figure('Name', 'ECGSYN Implicit Tustin');
plot(t, z, 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('z(t)');
title('ECGSYN - Fixed-step implicit Tustin');
grid on;

figure('Name', 'ECGSYN Implicit Tustin (zoom)');
idx_zoom = (t >= 5 & t <= 10);
plot(t(idx_zoom), z(idx_zoom), 'k', 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('z(t)');
title('ECGSYN - Fixed-step implicit Tustin (zoom)');
grid on;

% ============================================================
% Local functions
% ============================================================

function f = ecgsyn_rhs(x, tt, omega, ti, ai, bi, fhi, A)
    xx = x(1);
    yy = x(2);
    zz = x(3);

    r = sqrt(xx^2 + yy^2);
    alpha = 1 - r;

    theta = atan2(yy, xx);
    z0 = A * sin(2*pi*fhi*tt);

    zsum = 0.0;
    for k = 1:5
        dth = wrap_pm_pi(theta - ti(k));
        zsum = zsum + ai(k) * dth * exp(-(dth^2) / (2*bi(k)^2));
    end

    fx = alpha * xx - omega * yy;
    fy = alpha * yy + omega * xx;
    fz = -zsum - (zz - z0);

    f = [fx; fy; fz];
end

function J = ecgsyn_rhs_jacobian(x, tt, omega, ti, ai, bi, fhi, A)
    %#ok<INUSD>
    xx = x(1);
    yy = x(2);
    zz = x(3); %#ok<NASGU>

    r2 = xx^2 + yy^2;
    r  = sqrt(r2);

    if r < 1e-14
        r  = 1e-14;
        r2 = r^2;
    end

    % --------------------------------------------------------
    % x,y rows: exact analytic Jacobian
    % --------------------------------------------------------
    % f1 = x - x*r - omega*y
    % f2 = y - y*r + omega*x

    df1_dx = 1 - r - (xx^2)/r;
    df1_dy = -omega - (xx*yy)/r;
    df1_dz = 0;

    df2_dx =  omega - (xx*yy)/r;
    df2_dy = 1 - r - (yy^2)/r;
    df2_dz = 0;

    % --------------------------------------------------------
    % z row: exact analytic Jacobian
    % f3 = -sum_k ai * dth * exp(-dth^2/(2*bi^2)) - (z - z0)
    %
    % Let g(d) = d * exp(-d^2/(2b^2))
    % Then g'(d) = exp(-d^2/(2b^2)) * (1 - d^2/b^2)
    %
    % dth = wrap_pm_pi(theta - ti)
    % Away from branch cuts, d(dth)/d(theta) = 1
    % theta = atan2(y,x)
    % dtheta/dx = -y/(x^2+y^2)
    % dtheta/dy =  x/(x^2+y^2)
    % --------------------------------------------------------
    dtheta_dx = -yy / r2;
    dtheta_dy =  xx / r2;

    sum_dx = 0.0;
    sum_dy = 0.0;

    for k = 1:5
        dth = wrap_pm_pi(atan2(yy, xx) - ti(k));
        e   = exp(-(dth^2) / (2*bi(k)^2));
        gp  = e * (1 - (dth^2)/(bi(k)^2));   % derivative of d*exp(-d^2/(2b^2))

        sum_dx = sum_dx + ai(k) * gp * dtheta_dx;
        sum_dy = sum_dy + ai(k) * gp * dtheta_dy;
    end

    % f3 = -sum(...) - z + z0
    df3_dx = -sum_dx;
    df3_dy = -sum_dy;
    df3_dz = -1;

    J = [df1_dx, df1_dy, df1_dz;
         df2_dx, df2_dy, df2_dz;
         df3_dx, df3_dy, df3_dz];
end

function y = wrap_pm_pi(x)
    y = atan2(sin(x), cos(x));
end
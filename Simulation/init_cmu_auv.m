%% CMU-AUV Simulation Initialization Script
clear; clc;

% ==========================================
% 1. Physical Properties & Geometry (Table 1)
% ==========================================
m      = 11.95;                     % Mass [kg]
g      = 9.81;                     % Gravity [m/s^2] ~ 9.80998
W      = 117.23;                   % Weight [N]  (Table 3; = m*g with g = 9.81)
nabla  = 0.012138;                  % Volume [m^3]
B      = 119.071;                  % Buoyancy [N] (Table 3) -> slightly positive buoyancy, matches plant/EKF
rho    = 1000;                      % Water density [kg/m^3]

% Center of Gravity (CG) and Center of Buoyancy (CB) wrt Origin (CO)
r_g = [0; 0; 0.0620];               % [x_g; y_g; z_g] [m]
r_b = [0.0050; 0; 0.0213];          % [x_b; y_b; z_b] [m]

% Rigid-Body Inertia Tensor at CO [kg*m^2]
Ixx = 0.209849; Iyy = 0.305129; Izz = 0.264812;
Ixy = 0;        Ixz = 0.004125; Iyz = 0;

I_o = [ Ixx, -Ixy, -Ixz;
       -Ixy,  Iyy, -Iyz;
       -Ixz, -Iyz,  Izz ];

% Skew-symmetric matrix of r_g
S_rg = [   0,   -r_g(3),  r_g(2);
         r_g(3),    0,   -r_g(1);
        -r_g(2),  r_g(1),    0   ];

% Rigid-Body Mass Matrix (6x6)
M_RB = [ m * eye(3),    -m * S_rg;
         m * S_rg,        I_o    ];

% ==========================================
% 2. Hydrodynamic Added Mass (Table 3)
% ==========================================
M_A = [  5.55889,   0,       -0.02354,  0,        0.19245,   0;
         0,        13.33090,  0,       -0.78335,  0,         3.04845;
        -0.02972,   0,       15.01354,  0,       -3.43365,   0;
         0,        -0.82329,  0,        0.62734,  0,        -0.18743;
         0.20756,   0,       -3.42374,  0,        0.96326,   0;
         0,         3.05326,  0,       -0.17690,  0,         1.04373 ];

% Total System Inertia Matrix
M_total = M_RB + M_A;
inv_M   = inv(M_total);

% ==========================================
% 3. Damping Coefficients (Tables 4 & 5)
% ==========================================
D_L = [  0,        0,       -5.39181,  0,       -0.02484,  0;
         0,        0,        0,        0.08869,  0,       -0.00479;
         0.27142,  0,        1.13161,  0,       -0.03201,  0;
         0,        0.02142,  0,       -0.01501,  0,       -0.00040;
        -0.09545,  0,       -1.37277,  0,       -0.00294,  0;
         0,       -0.12307,  0,        0.01937,  0,       -0.01100 ];

% Quadratic damping (Table 6 of the System Overview PDF; rows X,Y,Z,K,M,N;
% columns u|u|, v|v|, w|w|, p|p|, q|q|, r|r|). NOTE: the plant, EKF and
% controller blocks use only abs(diag(D_Q)): 55.55 (surge), 137.3 (sway),
% 125.8 (heave), 1.17, 0.90, 1.33. (The previous version of this matrix had
% columns shuffled in rows 1,3,5,6, giving surge 1.30 and heave 0.)
D_Q = [ 55.5498,    0,        1.30135,  0,       -1.17923,  0;
         0,       137.29600,   0,        4.04258,  0,        0.02168;
        -6.90530,   0,       125.83700,  0,        0.53355,  0;
         0,       -10.77595,   0,       -1.17444,  0,        0.02980;
         5.47058,   0,       -28.77614,  0,       -0.90443,  0;
         0,        31.17960,   0,        0.91566,  0,       -1.33297 ];

% ==========================================
% 4. Thruster Configuration Matrix (Table 2)
% ==========================================
% Thruster mapping: tau = TCM * u_thr (u_thr: 8x1 thrust vector)
TCM = [  0.70711,  0.70711, -0.70711, -0.70711,  0.00000,  0.00000,  0.00000,  0.00000;
         0.70711, -0.70711,  0.70711, -0.70711,  0.00000,  0.00000,  0.00000,  0.00000;
         0.00000,  0.00000,  0.00000,  0.00000,  1.00000,  1.00000,  1.00000,  1.00000;
        -0.03111,  0.03111, -0.03111,  0.03111, -0.12000,  0.12000, -0.12000,  0.12000;
         0.03111,  0.03111, -0.03111, -0.03111, -0.05000, -0.05000,  0.05000,  0.05000;
         0.18385, -0.18385, -0.18385,  0.18385,  0.00000,  0.00000,  0.00000,  0.00000 ];

% Moore-Penrose Pseudo-Inverse for Allocation
TCM_pinv = pinv(TCM);

% Thruster limits (adjust according to thruster hardware specs)
T_max = 2.7*9.81;  % [N]
T_min = -2.4*9.81; % [N]

% ==========================================
% 5. SDRE / LQI Tuning Weights
% ==========================================
% State error vector: [e_x, e_y, e_z, e_phi, e_theta, e_psi, u, v, w, p, q, r, int_e_z]
Q_diag = [ 20, 20, 50, ...         % Position weights
           500,  500,  200, ...         % Orientation weights
           20,  20,   40, ...         % Linear velocity damping
           50,  50,   50, ...         % Angular velocity damping
           10 ];                     % Depth integral weight (LQI)
Q_base = diag(Q_diag);   % original single-set weights (kept for reference only)

% Control effort penalty (6-DOF generalized force/torque)
R_diag = [ 10, 10, 10, 50, 50, 50 ];
R_base = diag(R_diag);   % NOT used by the controller any more - see Q_tab / R_tab below


% ------------------------------------------
% Per-case LQR weights (used by the SDRE Controller block)
% ------------------------------------------
% Column k = diagonal of Q (13 states) / R (6 generalized forces) for
% trajectory case k:  1 = station keeping, 2 = straight line, 3 = circle,
% 4 = 3D infinity.  The SDRE Controller block picks the column from TrajCase,
% so changing TrajCase is all that is needed.  Only the ratio Q/R matters.
% Tuned with CMA-ES on a 20 Hz EKF + SDRE + feedforward simulation across
% three plant variants (nominal, +30% drag/+10% mass/+1.5 N buoyancy, and the
% opposite) so the gains are not fitted to a perfect model.
% Column 4 (infinity) reuses the circle gains: a full tuning run for it was
% not completed, but they beat the original weights on that case.
Q_tab = [
        80.04         1949         1359         1359    % x
        80.04         1949         1359         1359    % y
        6.966        9.682       0.6224       0.6224    % z
        674.7    2.043e+04         2977         2977    % phi
        674.7    2.043e+04         2977         2977    % theta
        267.6         2142         1397         1397    % psi
        535.2        18.15        7.989        7.989    % u
        535.2        18.15        7.989        7.989    % v
        440.5         15.3        603.8        603.8    % w
        4.406          654        1.429        1.429    % p
        4.406          654        1.429        1.429    % q
        18.15        12.93        113.5        113.5    % r
        986.7        18.52        16.19        16.19    % int_z
];

R_tab = [
        9.443       0.3461       0.5037       0.5037    % X
        9.443       0.3461       0.5037       0.5037    % Y
       0.6159       0.3895       0.3185       0.3185    % Z
        16.92         1007        726.5        726.5    % K
        16.92         1007        726.5        726.5    % M
        7.662        13.44        2.278        2.278    % N
];

% ==========================================
% 6. EKF Sensor Noise & Filter Covariance
% ==========================================
dt = 0.05; % 20 Hz control loop
Q_ekf = diag([ones(1,6)*1e-4, ones(1,6)*1e-2]);
% Measurement noise (1-sigma) used BOTH by the simulated sensors (Band-Limited
% White Noise blocks: Cov = sigma^2*dt) and by the EKF (R_ekf), so they match.
% These are assumed values - replace with your sensors' measured noise.
sigma_eta  = [0.03 0.03 0.01  0.01 0.01 0.01];  % [m m m rad rad rad] x,y (ultrasonic), z (depth), phi theta psi (IMU)
sigma_gyro = [0.01 0.01 0.01];                  % [rad/s] p q r (IMU gyro)
R_ekf = diag([sigma_eta.^2, sigma_gyro.^2]);

% ==========================================
% 7. Trajectory Case Selection
% ==========================================
%   1 = Station keeping
%   2 = Straight-line trajectory (constant depth)
%   3 = Circular trajectory (x-y plane, constant depth)
%   4 = 3D arbitrary infinity-shape
% TrajCase selects BOTH the reference (cmuauv_trajectory.m) and the tuned
% Q/R weights (column TrajCase of Q_tab / R_tab). Trajectory parameters
% (radius, speed, amplitudes, ...) are at the top of cmuauv_trajectory.m.
TrajCase = 4;

function x_ref = cmuauv_trajectory(t, traj_case)
%CMUAUV_TRAJECTORY  Reference generator for the CMU-AUV SDRE/LQI + EKF model.
%
%   x_ref = cmuauv_trajectory(t, traj_case)
%
%   traj_case : 1 = Station keeping
%               2 = Straight-line trajectory (constant depth)
%               3 = Circular trajectory (x-y plane, constant depth)
%               4 = 3D arbitrary infinity-shape (figure-eight with depth change)
%
%   x_ref = [x y z phi theta psi  u v w p q r  int_z]'   (13x1)
%     - position / Euler angles are in the NED frame (SNAME), z positive down
%     - u v w p q r are BODY-frame velocity references (consistent with the
%       controller error  x_est - x_ref, where x_est holds body velocities)
%     - int_z reference is 0
%
%   All tunable numbers for each case are in the parameter block below.
%   The vehicle is assumed to start at rest at the origin with zero attitude.
%   Yaw references are NOT wrapped: the plant yaw is an integral of r and is
%   unbounded, so a wrapped reference would create a 2*pi error step.

% ------------------------------------------------------------------ common
t_dive = 10.0;     % [s] time to reach the working depth / start pose (cases 1,2,4)

% ------------------------------------------------- case 1: station keeping
sk_x   = 0.0;      % [m]   target x (NED)
sk_y   = 0.0;      % [m]   target y (NED)
sk_z   = 1.0;      % [m]   target depth (positive down)
sk_psi = 0.0;      % [rad] target heading

% ----------------------------------------------- case 2: straight line
ln_depth   = 1.0;  % [m]   constant depth
ln_head    = 0.0;  % [rad] line heading in NED (0 = +x / North). Yaw follows it.
ln_length  = 5.0; % [m]   line length (starts at origin)
ln_speed   = 0.5;  % [m/s] cruise speed
ln_accel   = 0.1;  % [m/s^2] trapezoidal-profile acceleration
ln_tstart  = 10.0; % [s]   time the run along the line begins (>= t_dive)

% ----------------------------------------------- case 3: circle
%   Starts at the origin heading along +x at rest; after the dive the circle
%   phase rate is ramped 0 -> c_omega so there is no velocity step at the start.
c_R      = 1.0;    % [m]   radius (circle centre at (0, R))
c_omega  = 0.5;    % [rad/s] angular rate (one lap = 2*pi/omega = 12.6 s)
c_depth  = 1.5;    % [m]   constant depth
c_tramp  = 5.0;    % [s]   time to ramp the phase rate 0 -> omega
c_tstart = 10.0;   % [s]   time the circle begins (>= t_dive)

% ------------------------------------------ case 4: 3D infinity (figure-8)
%   x = inf_Ax*sin(th),  y = inf_Ay*sin(th)*cos(th),  z = inf_z0 - inf_Az*cos(th)
%   th = inf_omega*tau (smoothly accelerated over inf_tramp). Heading = path tangent.
inf_Ax    = 2.0;   % [m]   half-width  in x
inf_Ay    = 1.0;   % [m]   half-width  in y
inf_z0    = 1.5;   % [m]   mean depth
inf_Az    = 0.5;   % [m]   depth amplitude (path spans z0-Az ... z0+Az)
inf_omega = 0.3;   % [rad/s] phase rate (period 2*pi/omega = ~20.9 s)
inf_tramp = 5.0;   % [s]   time to ramp the phase rate 0 -> omega
inf_tstart = 10.0; % [s]   time the figure-eight begins (>= t_dive)

% ------------------------------------------------------------------ select
x_ref = zeros(13,1);
p   = zeros(3,1);   % NED position
pd  = zeros(3,1);   % NED velocity
psi = 0;  r = 0;    % heading, yaw rate

if traj_case == 3
    % ===== CASE 3: circle in the x-y plane at constant depth =====
    [h, dh] = traj_quintic(t, t_dive);              % dive to c_depth
    [th, thd] = traj_phase(t - c_tstart, c_omega, c_tramp);
    p   = [c_R*sin(th); c_R - c_R*cos(th); c_depth*h];
    pd  = [c_R*cos(th)*thd; c_R*sin(th)*thd; c_depth*dh];
    psi = th;                                       % tangent heading, unwrapped on purpose
    r   = thd;

elseif traj_case == 2
    % ===== CASE 2: straight line, trapezoidal speed profile =====
    [h, dh] = traj_quintic(t, t_dive);          % dive + align heading
    [s, sd] = traj_trapezoid(t - ln_tstart, ln_length, ln_speed, ln_accel);
    p   = [s*cos(ln_head); s*sin(ln_head); ln_depth*h];
    pd  = [sd*cos(ln_head); sd*sin(ln_head); ln_depth*dh];
    psi = ln_head*h;
    r   = ln_head*dh;

elseif traj_case == 4
    % ===== CASE 4: 3D infinity-shape =====
    tau = t - inf_tstart;
    [th, thd] = traj_phase(tau, inf_omega, inf_tramp);
    % heading of path tangent w.r.t. the phase (independent of speed)
    xp  = inf_Ax*cos(th);       yp  = inf_Ay*cos(2*th);     % d/dth
    xpp = -inf_Ax*sin(th);      ypp = -2*inf_Ay*sin(2*th);  % d2/dth2
    psi_raw = atan2(yp, xp);
    thm = mod(th, 2*pi);                    % unwrap (heading is periodic in th)
    if (thm > pi/2) && (thm < 3*pi/2) && (psi_raw > 0)
        psi_path = psi_raw - 2*pi;
    else
        psi_path = psi_raw;
    end
    dpsi_dth = (xp*ypp - yp*xpp)/(xp^2 + yp^2);
    % start pose (th = 0)
    psi0 = atan2(inf_Ay, inf_Ax);
    z0s  = inf_z0 - inf_Az;
    if tau <= 0
        [h, dh] = traj_quintic(t, t_dive);
        p   = [0; 0; z0s*h];
        pd  = [0; 0; z0s*dh];
        psi = psi0*h;
        r   = psi0*dh;
    else
        p   = [inf_Ax*sin(th); inf_Ay*sin(th)*cos(th); inf_z0 - inf_Az*cos(th)];
        pd  = [inf_Ax*cos(th)*thd; inf_Ay*cos(2*th)*thd; inf_Az*sin(th)*thd];
        psi = psi_path;
        r   = dpsi_dth*thd;
    end

else
    % ===== CASE 1 (default): station keeping =====
    [h, dh] = traj_quintic(t, t_dive);
    p   = [sk_x; sk_y; sk_z]*h;
    pd  = [sk_x; sk_y; sk_z]*dh;
    psi = sk_psi*h;
    r   = sk_psi*dh;
end

% NED velocity -> body-frame linear velocity (phi_d = theta_d = 0)
cps = cos(psi); sps = sin(psi);
u =  cps*pd(1) + sps*pd(2);
v = -sps*pd(1) + cps*pd(2);
w =  pd(3);

x_ref = [p(1); p(2); p(3); 0; 0; psi; u; v; w; 0; 0; r; 0];
end

% ======================================================================
function [h, dh] = traj_quintic(t, T)
% Minimum-jerk blend 0 -> 1 over [0,T]; dh is d(h)/dt.
s = min(max(t/T, 0), 1);
h = 10*s^3 - 15*s^4 + 6*s^5;
if t > 0 && t < T
    dh = (30*s^2 - 60*s^3 + 30*s^4)/T;
else
    dh = 0;
end
end

function [s, sd] = traj_trapezoid(tau, L, V, a)
% Distance s(tau) and speed sd(tau) of a trapezoidal profile covering L.
if L >= V^2/a
    Vp = V;  ta = V/a;  tc = (L - V^2/a)/V;
else
    Vp = sqrt(L*a);  ta = Vp/a;  tc = 0;      % triangular profile
end
T = 2*ta + tc;
if tau <= 0
    s = 0;  sd = 0;
elseif tau < ta
    s = 0.5*a*tau^2;  sd = a*tau;
elseif tau < ta + tc
    s = 0.5*a*ta^2 + Vp*(tau - ta);  sd = Vp;
elseif tau < T
    td = T - tau;
    s = L - 0.5*a*td^2;  sd = a*td;
else
    s = L;  sd = 0;
end
end

function [th, thd] = traj_phase(tau, omega, Tr)
% Phase th(tau) whose rate ramps linearly 0 -> omega over Tr, then is constant.
if tau <= 0
    th = 0;  thd = 0;
elseif tau < Tr
    th = omega*tau^2/(2*Tr);  thd = omega*tau/Tr;
else
    th = omega*(tau - Tr/2);  thd = omega;
end
end

%% ATTITUDE COMPUTATION
% Author: Laura Train
% Date 12/12
%% 
% Kalman filter to estimate position of a particle around a circle in the
% first quadrant
% Tool: NaveGo

matlabrc

addpath ../../kalman/

R = 2;
t = linspace(0,10,300)';

w = 2;


std_w = 0.001;
std_a = [0.002, 0.002];
std_pos = [0.01, 0.01];
std_v = [0.01, 0.01];

a_b = [-w^2*R, 0];

xr_noisy = R*ones(length(t),1) + std_pos(1)*randn(length(t),1);
xtheta_noisy = zeros(length(t),1) +  std_pos(2)*randn(length(t),1);
vr_noisy = zeros(length(t),1) +  std_v(1)*randn(length(t),1);
vtheta_noisy = w*R*ones(length(t),1) + std_v(2)*randn(length(t),1);

ar_noisy = -w^2*R*ones(length(t),1) +  std_a(1)*randn(length(t),1);
atheta_noisy = zeros(length(t),1) + std_a(2)*randn(length(t),1);


kf.xi = [R, 0, 0, w*R]';
kf.Pi = diag([0.05, 0.05, 0.05, 0.05]);
kf.z = [xr_noisy(1), xtheta_noisy(1), vr_noisy(1), vtheta_noisy(1)]';

syms x y vx vy 

DCMbn(1,1) = cos(atan(y/x));
DCMbn(1,2) = -sin(atan(y/x));
DCMbn(2,1) = sin(atan(y/x));
DCMbn(2,2) = cos(atan(y/x));

DCMnb(1,1) = cos(atan(y/x));
DCMnb(1,2) = sin(atan(y/x));
DCMnb(2,1) = -sin(atan(y/x));
DCMnb(2,2) = cos(atan(y/x));

O = zeros(2);
h = [DCMnb O; O DCMnb]*[x; y; vx; vy];

H = jacobian(h, [x, y, vx, vy]);

kf.H = double(subs(H, [x, y, vx, vy], ...
                      [kf.xi(1), kf.xi(2), kf.xi(3), kf.xi(4)]));

kf.R = diag([std_pos(1), std_pos(2), std_v(1), std_v(2)]).^2;
kf = kf_update(kf);
kf.G = eye(4);
kf.Q = diag([std_pos(1), std_pos(2), std_v(1), std_v(2)]).^2;
state(1,1:4) = kf.xp;


for i = 2:27
    dt = t(i) - t(i-1);
    kf.H = double(subs(H, [x, y, vx, vy], ...
                      [kf.xi(1), kf.xi(2), kf.xi(3), kf.xi(4)]));

    kf.z = [xr_noisy(i), xtheta_noisy(i), vr_noisy(i), vtheta_noisy(i)]';
    kf.F = [1, 0, dt, 0;
            0, 1, 0, dt;
            0, 0, 1, 0;
            0, 0, 0, 1];
    kf.B = [1/2*dt^2, 0;
            0, 1/2*dt^2;
            dt, 0;
            0, dt];
    kf.u = [cos(atan2(kf.xp(2),kf.xp(1))), -sin(atan2(kf.xp(2),kf.xp(1)));
            sin(atan2(kf.xp(2),kf.xp(1))), cos(atan2(kf.xp(2),kf.xp(1)))]*[ar_noisy(i); atheta_noisy(i)];
        
    %kf_prediction(kf, dt)
    kf = kalman_circ(kf);
    state(i,:) = [kf.xp(1), kf.xp(2), kf.xp(3), kf.xp(4)];
    
end

x = state(:,1);
y = state(:,2);
plot(x,y, '*r')
axis equal
xlabel('x')
ylabel('y')
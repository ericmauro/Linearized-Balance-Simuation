% Eric Mauro
% Balance Control Dynamic Model
%   Model derived from 3-segment inverted pendulum w/ human constants,
%   linearized around about upright position
%   joint angles converted to segment angles
% Created: Jan-27-2017
% Updated: Mar-31-2017

function []=balance_dynamics5()
clc; clear all; close all;

noise = 0;

%% Physical Parameters
g = 9.81; % Gravitational constant (m/s/s)

m = [4;          % Mass of calf/lower leg (kg)
    7;           % Mass of thigh (kg)
    48];         % Mass of torso + head (kg)

L = [0.6;        % Length of calf/lower leg (m)
    0.5;         % Length of thigh (m)
    0.8];        % Length of torso + head (m)

c = [0.567;      % Ratio of l1, distance from ankle to CoM for calf
    0.567;       % Ratio of l2, distance from knee to CoM for thigh
    0.4829];     % Ratio of l3, distance to CoM of torso + head

I = [0.12;       % Moment of Inertia for calf (kg*m^2)
    0.1458;      % Moment of Inertia for thigh (kg*m^2)
    2.26];       % Moment of Inertia for torso + head (kg*m^2)

% Range of motion angle constraints (ankle, knee, hip)
rom = deg2rad([-50 20; -135 0; -30 120]);

%% Linearized Model
% Linearized around upright position (theta = 0)
% cos -> 1, sin -> theta, \dot{theta}^2 -> 0
h = [m(1)*c(1)^2+m(2)+m(3); m(2)*c(2)^2+m(3); m(3)*c(3)^2]; % Constant
k = [m(1)*c(1)+m(2)+m(3); m(2)*c(2)+m(3); m(3)*c(3)];       % Constant

M = [I(1)+h(1)*L(1)^2 L(1)*L(2)*k(2) L(1)*L(3)*k(3);
    L(1)*L(2)*k(2) I(2)+h(2)*L(2)^2 L(2)*L(3)*k(3);
    L(1)*L(3)*k(3) L(2)*L(3)*k(3) I(3)+h(3)*L(3)^2]; % \ddot{theta} coeff. matrix

G = g.*[L(1)*k(1) 0 0; 
        0 L(2)*k(2) 0; 
        0 0 L(3)*k(3)]; % theta coeff. matrix
    
% Convert to segment angles
S = [1 0 0; 1 1 0; 0 1 1]; % Transformation Matrix
%Ms = [M(1,1)+M(1,2) M(1,2)+M(1,3) M(1,3);
      %M(1,2)+M(2,2) M(2,2)+M(2,3) M(2,3);
      %M(1,3)+M(2,3) M(2,3)+M(3,3) M(3,3)];
Ms = M*S;
  
%Gs = [G(1,1) 0 0; G(2,2) G(2,2) 0; 0 G(3,3) G(3,3)];
Gs = G*S;

%% Set up State-Space Model and Simulation
A = [zeros(3) eye(3); inv(Ms)*Gs zeros(3)];
B = [zeros(3); inv(Ms)];
C = [eye(3) zeros(3,3)];
n = 6;

R = eye(3);
%Q1 = 100*diag([1,5,20]);
%Q2 = 10*diag([30,5,20]);
Q1 = diag([4900,18225,22500]);
Q2 = 0.25*diag([40000,18225,22500]);
Q = [Q1 zeros(3); zeros(3) Q2];

[Ko,Po,Eo] = lqr(A,B,Q,R);

x0 = [0.1 -0.3 0.1 0 0 0]; % Lean Back
%x0 = [0.25 -1.45 1.75 0 0 0]; % Squat

%% Simulate and Plot
dt = 0.02;          % Time step
T = 2;              % Time period
N = T/dt;           % Number of time steps
X = [x0,x0];        % Initialize both models
x_save = [];
t_save = [];
u_save = [];

for i = 1:N
    [t,X] = ode15s(@pendmodel,[(i-1)*dt i*dt],X(end,:));
    x_save=[x_save;X];
    t_save=[t_save;t];
end

% Generate exploratory signal to plot input signals
if noise == 1
    w1 = 1;
    w2 = 10;
    w3 = 100;
else
    w1 = 0;
    w2 = 0;
    w3 = 0;
end
e = 50.*[(sin(w1*t_save')); (sin(w2*t_save')); (sin(w3*t_save'))]; % Exploratory signal
u_n = -Ko*x_save(:,1:n)'+e;
u_l = -Ko*x_save(:,n+1:end)'+e;

% Plot State Data
figure
subplot(3,1,1)
hold on
plot(t_save,x_save(:,1))
plot(t_save,zeros(length(t_save),1),'r--')
plot(t_save,x_save(:,n+1),'g')
hold off
xlabel('Time (s)')
ylabel('\phi_1')

subplot(3,1,2)
hold on
plot(t_save,x_save(:,2))
plot(t_save,zeros(length(t_save),1),'r--')
plot(t_save,x_save(:,n+2),'g')
hold off
xlabel('Time (s)')
ylabel('\phi_2')

subplot(3,1,3)
hold on
plot(t_save,x_save(:,3))
plot(t_save,zeros(length(t_save),1),'r--')
plot(t_save,x_save(:,n+3),'g')
hold off
xlabel('Time (s)')
ylabel('\phi_3')

% Plot input data
figure
subplot(3,1,1)
plot(t_save,u_n)
xlabel('Time (s)')
ylabel('Input, u (Nonlinear model)')
subplot(3,1,2)
plot(t_save,u_l)
xlabel('Time (s)')
ylabel('Input, u (Linear model)')
subplot(3,1,3)
plot(t_save,abs(u_n-u_l))
xlabel('Time (s)')
ylabel('|u_n - u_l|')

%% Pendulum Simulation
x1 = zeros(length(t_save),1); % Ankle Stationary
y1 = zeros(length(t_save),1);
x2 = x1+L(1)*sin(x_save(:,1)); % Knee position
y2 = y1+L(1)*cos(x_save(:,1));
x3 = x2+L(2)*sin(x_save(:,1)+x_save(:,2)); % Hip Position
y3 = y2+L(2)*cos(x_save(:,1)+x_save(:,2));
x4 = x3+L(3)*sin(x_save(:,2)+x_save(:,3)); % Head Position
y4 = y3+L(3)*cos(x_save(:,2)+x_save(:,3));

figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

for i = 1:50:length(t_save)
   subplot(3,3,[1:6])
   h=plot([x1(i) x2(i) x3(i) x4(i)],[y1(i) y2(i) y3(i) y4(i)],'MarkerSize',30,'Marker','.','LineWidth',2);
   axis([-2 2 0 3])
   xlabel('Horizontal Distance (m)')
   ylabel('Vertical Distance (m)')
   
   % Add Face
   axes('pos',[0.25*x4(i)+0.325 y4(i)/2.67 .4 .2])
   [im, map, alpha] = imread('JGLFace13.png');
   f = imshow(im);
   set(f, 'AlphaData', alpha);
   

   F(i) = getframe;

   subplot(3,3,7)
   hold on
   plot(t_save(1:i),x_save(1:i,1))
   plot(t_save(1:i),zeros(i,1),'r--')
   hold off
   xlabel('Time (s)')
   ylabel('\phi_1')

   subplot(3,3,8)
   hold on
   plot(t_save(1:i),x_save(1:i,2))
   plot(t_save(1:i),zeros(i,1),'r--')
   hold off
   xlabel('Time (s)')
   ylabel('\phi_2')

   subplot(3,3,9)
   hold on
   plot(t_save(1:i),x_save(1:i,3))
   plot(t_save(1:i),zeros(i,1),'r--')
   hold off
   xlabel('Time (s)')
   ylabel('\phi_3')
end

function dX = pendmodel(t,X)
    x = X(1:n);
    xl = X(n+1:end);
    % Nonlinear model with linear feedback control (LQR) on input
    Mn = [I(1)+h(1)*L(1)^2 L(1)*L(2)*k(2)*cos(-x(2)) L(1)*L(3)*k(3)*cos(x(1)-x(2)-x(3));
          L(1)*L(2)*k(2)*cos(-x(2)) I(2)+h(2)*L(2)^2 L(2)*L(3)*k(3)*cos(x(1)-x(3));
          L(1)*L(3)*k(3)*cos(x(1)-x(2)-x(3)) L(2)*L(3)*k(3)*cos(x(1)-x(3)) I(3)+h(3)*L(3)^2]; 

    Nn = [0 L(1)*L(2)*k(2)*sin(-x(2)) L(1)*L(3)*k(3)*sin(x(1)-x(2)-x(3));
          -L(1)*L(2)*k(2)*sin(-x(2)) 0 L(2)*L(3)*k(3)*sin(x(1)-x(3));
          -L(1)*L(3)*k(3)*sin(x(1)-x(2)-x(3)) -L(2)*L(3)*k(3)*sin(x(1)-x(3)) 0];

    Gn = g.*[L(1)*k(1)*sin(x(1)); 
             L(2)*k(2)*sin(x(2)+x(1)); 
             L(3)*k(3)*sin(x(3)+x(2))];
    
    % Adjust M and N for segment angle coordinate system
    Mns = [Mn(1,1)+Mn(1,2) Mn(1,2)+Mn(1,3) Mn(1,3);
           Mn(1,2)+Mn(2,2) Mn(2,2)+Mn(2,3) Mn(2,3);
           Mn(1,3)+Mn(2,3) Mn(2,3)+Mn(3,3) Mn(3,3)];
    
    Nns1 = [Nn(1,2) Nn(1,2)+Nn(1,3) Nn(1,3);
            -Nn(1,2) Nn(2,3) Nn(2,3);
            -(Nn(1,3)+Nn(2,3)) -Nn(2,3) 0];
            
    Nns2 = 2*x(5).*[Nn(1,2)*x(4)+Nn(1,3)*x(6);
                    Nn(2,3)*x(6);
                    Nn(2,3)*x(4)];
    
    % Add exploratory noise to input
    if noise == 1
        w1 = 1;
        w2 = 10;
        w3 = 100;
    else
        w1 = 0;
        w2 = 0;
        w3 = 0;
    end
    e = 50.*[(sin(w1*t)); (sin(w2*t)); (sin(w3*t))]; % Exploratory signal
    u = -Ko*x+e;
    
    X1 = inv(Mns)*Nns1*(x(4:6).^2);
    X2 = inv(Mns)*Nns2;
    X3 = inv(Mns)*Gn;
    X4 = inv(Mns)*u;

    dx = [x(4:6); X1+X2+X3+X4];
    
    % Linearized
    ul = -Ko*xl+e;
    dxl = A*xl + B*ul;
    
    dX = [dx;dxl];
end
end

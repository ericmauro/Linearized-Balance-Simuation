% Eric Mauro
% Balance Control Dynamic Model
%   Model derived from 3-segment inverted pendulum w/ human constants,
%   linearized around about upright position
% Created: Jan-27-2017
% Updated: Feb-16-2017
clc; clear all; close all;
control = 1; % Open-loop (0) or Feedback (1)

%% Physical Parameters
g = 9.81; % Gravitational constant (m/s/s)

m = [3.52;          % Mass of calf/lower leg (kg)
    7.77;           % Mass of thigh (kg)
    34.66+6.07];    % Mass of torso + head (kg)

l = [0.43;          % Length of calf/lower leg (m)
    0.40;           % Length of thigh (m)
    0.52+0.19];     % Length of torso + head (m)

c = [0.567;         % Ratio of l1, distance from ankle to CoM for calf
    0.567;          % Ratio of l2, distance from knee to CoM for thigh
    0.4829];        % Ratio of l3, distance to CoM of torso + head

I = [0.0088;        % Moment of Inertia for calf (kg*m^2)
    0.0200;         % Moment of Inertia for thigh (kg*m^2)
    2.17+0.02];     % Moment of Inertia for torso + head (kg*m^2)

% Range of motion angle constraints (ankle, knee, hip)
rom = deg2rad([-50 20; -135 0; -30 120]);

%% Linearized Model
% Linearized around upright position (theta = 0)
% cos -> 1, sin -> theta, \dot{theta}^2 -> 0
h = [m(1)*c(1)^2+m(2)+m(3); m(2)*c(2)^2+m(3); m(3)*c(3)^2]; % Constant
k = [m(1)*c(1)+m(2)+m(3); m(2)*c(2)+m(3); m(3)*c(3)];       % Constant

M = [I(1)+h(1)*l(1)^2 l(1)*l(2)*k(2) l(1)*l(3)*k(3);
    l(1)*l(2)*k(2) I(2)+h(2)*l(2)^2 l(2)*l(3)*k(3);
    l(1)*l(3)*k(3) l(2)*l(3)*k(3) I(3)+h(3)*l(3)^2]; % \ddot{theta} coeff. matrix

P = g.*[l(1)*k(1) 0 0; 
        0 l(2)*k(2) 0; 
        0 0 l(3)*k(3)]; % theta coeff. matrix

%% Set up State-Space Model and Simulation
A = [zeros(3) eye(3); inv(M)*P zeros(3)];
B = [zeros(3); inv(M)];
C = [eye(3) zeros(3,3)];

R = eye(3);
%Q1 = eye(3);
Q1 = 10000*diag([1,5,20]);
Q2 = 1000*diag([30,5,20]);
Q = [Q1 zeros(3); zeros(3) Q2];

[Ko,Po,Eo] = lqr(A,B,Q,R);

if control == 1 % Feedback control simulation
    sys = ss(A-B*Ko,B,C,0);
elseif control == 0 % Open-loop simulation
    sys = ss(A,B,C,0);
end

t = 0:0.02:2; % Time
ut = zeros(3,length(t)); % Zero-input/Open-loop
x0 = [0.25 -1.2 0.3 0 0 0]; % Initial states
%x0 = [-0.4 -0.45 -0.45 0 0 0];
%x0 = [0.2 -0.2 1.5 0 0 0];

%% Simulate and Plot
Y = lsim(sys,ut,t,x0);

% Saturation
for i = 1:3
    Y(:,i) = Y(:,i).*((rom(i,1)<Y(:,i))&(Y(:,i)<rom(i,2))) + rom(i,1).*(Y(:,i)<=rom(i,1)) + rom(i,2).*(Y(:,i)>=rom(i,2));
end

x1 = zeros(length(t),1); % Ankle Stationary
y1 = zeros(length(t),1);
x2 = x1+l(1)*sin(Y(:,1)); % Knee position
y2 = y1+l(1)*cos(Y(:,1));
x3 = x2+l(2)*sin(Y(:,2)); % Hip Position
y3 = y2+l(2)*cos(Y(:,2));
x4 = x3+l(3)*sin(Y(:,3)); % Head Position
y4 = y3+l(3)*cos(Y(:,3));

figure;
pause(0.00001);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);

for i = 1:length(t)
   subplot(3,3,[1:6])
   h=plot([x1(i) x2(i) x3(i) x4(i)],[y1(i) y2(i) y3(i) y4(i)],'MarkerSize',30,'Marker','.','LineWidth',2);
   axis([-2 2 0 2])
   xlabel('Horizontal Distance (m)')
   ylabel('Vertical Distance (m)')
   %{
   % Add Face
   axes('pos',[0.25*x4(i)+0.325 y4(i)/2 .4 .2])
   [im, map, alpha] = imread('JGLFace13.png');
   f = imshow(im);
   set(f, 'AlphaData', alpha);
   %}

   F(i) = getframe;

   subplot(3,3,7)
   hold on
   plot(t(1:i),Y(1:i,1))
   plot(t(1:i),zeros(i,1),'r--')
   hold off
   xlabel('Time (s)')
   ylabel('\theta_1')

   subplot(3,3,8)
   hold on
   plot(t(1:i),Y(1:i,2))
   plot(t(1:i),zeros(i,1),'r--')
   hold off
   xlabel('Time (s)')
   ylabel('\theta_2')

   subplot(3,3,9)
   hold on
   plot(t(1:i),Y(1:i,3))
   plot(t(1:i),zeros(i,1),'r--')
   hold off
   xlabel('Time (s)')
   ylabel('\theta_3')
end

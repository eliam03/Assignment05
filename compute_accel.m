%Computes the linear and angular acceleration of the box
%given its current position and orientation

%INPUTS:
%x: current x position of the box
%y: current y position of the box
%theta: current orientation of the box

%box_params: a struct containing the parameters that describe the system
%Fields:
%box_params.m: mass of the box
%box_params.I: moment of inertia w/respect to centroid
%box_params.g: acceleration due to gravity
%box_params.k_list: list of spring stiffnesses
%box_params.l0_list: list of spring natural lengths
%box_params.P_world: 2 x n list of static mounting
% points for the spring (in the world frame)
%box_params.P_box: 2 x n list of mounting points
% for the spring (in the box frame)


%OUTPUTS
%ax: x acceleration of the box
%ay: y acceleration of the box
%atheta: angular acceleration of the box
function [ax,ay,atheta] = compute_accel(x,y,theta,box_params)

%Deconstruct the struct
m = box_params.m; % mass of the box
I = box_params.I; % moment of inertia w/respect to centroid
g = box_params.g; % acceleration due to gravity
k = box_params.k_list; % list of spring stiffnesses
l0 = box_params.l0_list; % list of spring natural lengths
P_world = box_params.P_world; % 2 x n list of static mounting points for the spring (in the world frame)
P_box = box_params.P_box; % 2 x n list of mounting points for the spring (in the box frame)

%Center of mass
P_c = [x ; y];
% Calculate the Spring Force (might be wrong)
F_spring = compute_spring_force(k,l0,P_box(1,:),P_box(2,:));

% Calculate the external force acting on the system
F_ext = [0; m*g] + F_spring;

% Calculate the external torque (might be wrong)
tau_ext = cross(P_world - P_c, F_spring);

% Calulate accelerations
a = F_ext./m;
ax = a(1);
ay = a(2);
atheta = tau_ext./I;
end

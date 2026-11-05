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
    
    % Deconstruct the struct
    m = box_params.m; % mass of the box
    I = box_params.I; % moment of inertia w/respect to centroid
    g = box_params.g; % acceleration due to gravity
    k_list = box_params.k_list; % list of spring stiffnesses
    l0_list = box_params.l0_list; % list of spring natural lengths
    P_world = box_params.P_world; % 2 x n list of static mounting points for the spring (in the world frame)
    P_box = box_params.P_box; % 2 x n list of mounting points for the spring (in the box frame)
    
    % Center of mass
    P_c = [x ; y];
    
    % Transform box points into world coordinates 
    P_box_world = compute_rbt(x, y, theta, P_box);

    % Initial force and torque 
    F_grav = [0; m*g]; % Gravity
    tau_ext = 0; % No initial external torque
    
    % Simulate for every spring

    for i = 1:length(k_list)
        % Fixed end of the spring
        PA = P_world(:,i);

        % Moving end of the spring 
        PB = P_box_world(:,i);

        % Calculate the Spring Force (might be wrong)
        F_i = compute_spring_force(k_list(i),l0_list(i),PA,PB);

        % Update the external force by adding the spring force
        F_ext = F_grav + F_i;

        % Position from the centroid to the spring force application point
        r_i = PB - P_c;

        % Torque exerted by each spring
         tau_i = r_i(1)*F_i(2) - r_i(2)*F_i(1);
        
        % Add the spring torque
        tau_ext = tau_ext + tau_i;
    end 

    % Calulate accelerations
    ax = F_ext(1)/m;
    ay = F_ext(2)/m;
    atheta = tau_ext/I;
end

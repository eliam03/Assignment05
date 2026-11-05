function simulate_box()
    %define system parameters
    box_params = struct();
    box_params.m = %your code here
    box_params.I = %your code here
    box_params.g = %your code here
    box_params.k_list = %your code here
    box_params.l0_list = %your code here
    box_params.P_world = %your code here
    box_params.P_box = %your code here
    %load the system parameters into the rate function
    %via an anonymous function
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    x0 = %your code here
    y0 = %your code here
    theta0 = %your code here
    vx0 = %your code here
    vy0 = %your code here
    vtheta0 = %your code here
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = %your code here
    %run the integration
        DormandPrince = struct();

DormandPrince.C = ...
    [0, 1/5, 3/10, 4/5, 8/9, 1, 1];

DormandPrince.B = ...
    [35/384, 0, 500/1113, 125/192, ...
     -2187/6784, 11/84, 0;
     5179/57600, 0, 7571/16695, 393/640, ...
     -92097/339200, 187/2100, 1/40];

DormandPrince.A = ...
    [0, 0, 0, 0, 0, 0, 0;
     1/5, 0, 0, 0, 0, 0, 0;
     3/40, 9/40, 0, 0, 0, 0, 0;
     44/45, -56/15, 32/9, 0, 0, 0, 0;
     19372/6561, -25360/2187, 64448/6561, ...
        -212/729, 0, 0, 0;
     9017/3168, -355/33, 46732/5247, ...
        49/176, -5103/18656, 0, 0;
     35/384, 0, 500/1113, 125/192, ...
        -2187/6784, 11/84, 0];
    % [tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...)
end
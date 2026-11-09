
function simulate_box_equalibrium()
    clf;
    close all;
    %define system parameters
    LW = 10; LH = 10; LG = 3;
    m = 20; Ic = (1/12)*(LH^2+LW^2);
    g = 9.1; k = 50; k_list = [.5*k,.5*k,2*k,5*k];
    l0 = 1.5*LG;
    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;
    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];
    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];
    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];
    P_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];
    P_box = [Ptl_box,Pbl_box,Ptr_box,Pbr_box];
    %define system parameters
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(P_world,2));
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;
    %load the system parameters into the rate function
    %via an anonymous function
   
    x0 = 5;
    y0 = 5;
    theta0 = 45;
    vx0 = 2;
    vy0 = 3;
    vtheta0 = 35;
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = [0, 100];
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

    t=3;
    solver_params=struct();

    box_wrapper = @(V0) box_rate_func(t,V0,box_params);

    h = 0.05;
    p = 5.5;
    error_desired = 0.0001;

    [Xn,~] = multi_newton_solver(box_wrapper,V0,solver_params);
    
    dx0 = Xn(1);
    dy0 = Xn(2);
    dtheta0 = Xn(3);
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;
    %small number used to scale initial perturbation
    epsilon =0.00005;

    Veq=Xn;

    V0 = Veq + epsilon*[dx0;dy0;dtheta0;vx0;vy0;vtheta0];

    my_rate_func = @(t_in, V_in) box_rate_func(t_in,V_in,box_params);

    equil_func = @(V_in) box_rate_func(0, V_in, box_params);
    j_approx = approximate_jacobian(equil_func, Veq)

    my_linear_rate = @(t_in, V_in) j_approx*(V_in-Veq);

    %run the integration of nonlinear system
    [tlist_nonlinear,Vlist_nonlinear] =explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h,DormandPrince,p,error_desired);
    %run the integration of linear system
    [tlist_linear,Vlist_linear] = explicit_RK_variable_step_integration(my_linear_rate,tspan,V0,h,DormandPrince,p,error_desired);
    
    subplot(3, 1, 1)
    plot(tlist_nonlinear,Vlist_nonlinear(:,1))
    plot(tlist_linear,Vlist_linear(:,1))
    % x
    subplot(3, 1, 2)
    plot(tlist_nonlinear,Vlist_nonlinear(:, 2))
    plot(tlist_nonlinear,Vlist_nonlinear(:, 2))
    % y
    subplot(3, 1, 3)
    plot(tlist_nonlinear,Vlist_nonlinear(:, 3))
    plot(tlist_nonlinear,Vlist_nonlinear(:, 3))
    % theta
end
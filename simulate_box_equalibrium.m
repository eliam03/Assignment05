
function simulate_box_equalibrium()
    clf;
    close all;

    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

    %define system parameters
  LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2+LW^2);
    g = 1; k = 20; k_list = [.5*k,.5*k,2*k,5*k];
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
    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];

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
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);

    x0 = .182;
    y0 = .478;
    theta0 = -.321;
    vx0 = 0.5;
    vy0 = 0.5;
    vtheta0 = -.5;
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = [0, 60];
    t=3;
    solver_params=struct();

    box_wrapper = @(Vin) box_rate_func(t,Vin,box_params);

    h = 0.05;
    p = 5.5;
    error_desired = 0.0001;

    [Xn,~] = multi_newton_solver(box_wrapper,V0,solver_params);
    
    box_wrapper(Xn)
    
    dx0 = Xn(1);
    dy0 = Xn(2);
    dtheta0 = Xn(3);
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;
    %small number used to scale initial perturbation
    epsilon =10;

    Veq=Xn;

    V0 = Veq + epsilon*[dx0;dy0;dtheta0;vx0;vy0;vtheta0];

    my_rate_func = @(t_in, V_in) box_rate_func(t_in,V_in,box_params);

    equil_func = @(V_in) box_rate_func(0, V_in, box_params);
    j_approx = approximate_jacobian(equil_func, Veq);

    my_linear_rate = @(t_in, V_in) j_approx*(V_in-Veq);

    %run the integration of nonlinear system
    % [tlist_nonlinear,Vlist_nonlinear] =explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h,DormandPrince,p,error_desired);
    [tlist_nonlinear,Vlist_nonlinear] = ode45(my_rate_func,tspan,V0);
    %run the integration of linear system
    % [tlist_linear,Vlist_linear] = explicit_RK_variable_step_integration(my_linear_rate,tspan,V0,h,DormandPrince,p,error_desired);
    [tlist_linear,Vlist_linear] = ode45(my_linear_rate,tspan,V0);
    
   figure();
   t = tiledlayout(3, 1);
    % subplot(3, 1, 1)
    ax1 = nexttile;
    plot(ax1, tlist_nonlinear,Vlist_nonlinear(:,1)); hold on;
    plot(ax1, tlist_linear,Vlist_linear(:,1),'r--')
    title("Large Perturbation of a System in Equilibrium Over Time ($\varepsilon$ = 10)")
    xlabel("Time (-)")
    ylabel("X Position (-)")
    hold off;
    % x
    % subplot(3, 1, 2)
    ax2 = nexttile;
    plot(ax2, tlist_nonlinear,Vlist_nonlinear(:, 2)); hold on;
    plot(ax2, tlist_linear,Vlist_linear(:, 2),"r--"); hold off;
    xlabel("Time (-)")
    ylabel("Y Position (-)")
    % y
    % subplot(3, 1, 3)
    ax3 = nexttile;
    plot(ax3, tlist_nonlinear,Vlist_nonlinear(:, 3)); hold on;
    plot(ax3, tlist_linear,Vlist_linear(:, 3),"r--"); hold off;
    xlabel("Time (-)")
    ylabel("Theta (-)")
    % theta
    legend("Non linear", "Linear", 'location', 'southoutside', 'Orientation','horizontal', Interpreter="latex")
    sgtitle('Comparison of linearized and nonlinear ODE ', Interpreter="latex")

    myfig = gcf();
    exportgraphics(myfig , 'large_perturbation.png', 'Resolution', 600);
end

function simulate_box()
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
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
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

    box_wrapper = @(t, V0) box_rate_func(t,V0,box_params);
    h = 0.05;
    p = 5.5;
    error_desired = 0.0001;

    [tlist,Vlist, ~, ~, ~] = explicit_RK_variable_step_integration(box_wrapper,tspan,V0,h,DormandPrince,p,error_desired);
    

    fname='box_with_springs_animation.avi';
    input_fname = fname;
    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    open(writerObj); %must call open before writing any frames
    fig1 = figure(1);
    xlabel("$X$ Position (-)", Interpreter="latex");
    ylabel("$Y$ Position (-)", Interpreter="latex");
    title("Box Simulation", Interpreter="latex")
    set(gca,'TickLabelInterpreter','latex')
    fontsize(scale=2)

    set(fig1,'units','pixels','position',[0 0 1440 1080])
    num_zigs = 5;
    w = .1;
    hold on;
    % lowk inspired by shauna's code.... 
    sp1 = initialize_spring_plot(num_zigs,w);
    sp2 = initialize_spring_plot(num_zigs,w);
    sp3 = initialize_spring_plot(num_zigs,w);
    sp4 = initialize_spring_plot(num_zigs,w);
    origin = plot(0, 0, 'black');
    axis equal; axis square;
    axis([-40,40,-40,40]);
    for frames = 1:200
        x = Vlist(frames,1);
        y = Vlist(frames, 2);
        theta = Vlist(frames, 3);
        P_box = box_params.P_box;
        P_world = compute_rbt(x, y, theta, P_box);
        Q_world = box_params.P_world;
        update_spring_plot(sp1,P_world(:, 1), Q_world(:, 1));
        update_spring_plot(sp2,P_world(:, 2), Q_world(:, 2));
        update_spring_plot(sp3,P_world(:, 3), Q_world(:, 3));
        update_spring_plot(sp4,P_world(:, 4), Q_world(:, 4));

        %Draw box
        bound = box_params.boundary_pts;
        box_bound = compute_rbt(x, y, theta, bound);
        set(origin, 'xdata', box_bound(1, :), 'ydata', box_bound(2, :));
        drawnow;

        %capture a frame (what is currently plotted)
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);

        
    end
    close(writerObj);
end
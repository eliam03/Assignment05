
function simulate_box()
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

    
    box_wrapper = @(t, V0) box_rate_func(t,V0,box_params);


    [tlist,Vlist] = ode45(box_wrapper,linspace(0,tspan(end),1001),V0);
    

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

    sp1 = initialize_spring_plot(num_zigs,w);
    sp2 = initialize_spring_plot(num_zigs,w);
    sp3 = initialize_spring_plot(num_zigs,w);
    sp4 = initialize_spring_plot(num_zigs,w);
    origin = plot(0, 0, 'black');
    axis equal; axis square;
    axis([-40,40,-40,40]);
    for frames = 1:length(tlist)
        x = Vlist(frames,1);
        y = Vlist(frames, 2);
        theta = Vlist(frames, 3);
        R_world = compute_rbt(x, y, theta, box_params.P_box);
        Q_world = box_params.P_world;
        update_spring_plot(sp1,R_world(:, 1), Q_world(:, 1));
        update_spring_plot(sp2,R_world(:, 2), Q_world(:, 2));
        update_spring_plot(sp3,R_world(:, 3), Q_world(:, 3));
        update_spring_plot(sp4,R_world(:, 4), Q_world(:, 4));

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
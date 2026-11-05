function spring_plotting_example()
    num_zigs = 5;
    w = .1;
    hold on;
    spring_plot_struct = initialize_spring_plot(num_zigs,w);
    axis equal; axis square;
    axis([-3,3,-3,3]);
    for theta=linspace(0,6*pi,1000)
        P1 = [.5;.5];
        P2 = 2*[cos(theta);sin(theta)];
        update_spring_plot(spring_plot_struct,P1,P2)
        drawnow;
    end
end
function [XB1,XB2,num_evals, diff] = RK_step_embedded( rate_func_in,t,XA,h,BT_struct)
    diff = norm(rate_func_in(t + h, XA) - rate_func_in(t, XA));
    % Extract Butcher tableau
    A = BT_struct.A;
    B = BT_struct.B;
    C = BT_struct.C(:);

    % Number of RK stages
    num_stages = size(A,1);

    % Store k1, k2, ..., ks
    K = zeros( ...
        numel(XA), ...
        num_stages);


    %% Calculate RK stages

    for i = 1:num_stages

        if i == 1

            sum_val = zeros(size(XA));

        else

            sum_val = ...
                K(:,1:i-1) ...
                * A(i,1:i-1)';

        end


        K(:,i) = ...
            rate_func_in( ...
            t + C(i)*h, ...
            XA + h*sum_val);

    end

    % Dormand-Prince 5th-order estimate
    XB1 = XA + h*(K*B(1,:)');


    % Dormand-Prince 4th-order estimate
    XB2 = XA + h*(K*B(2,:)');


    % One rate-function evaluation per stage
    num_evals = num_stages;

end
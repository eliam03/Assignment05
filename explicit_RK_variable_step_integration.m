%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, step_fail_rate] = explicit_RK_variable_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)
    num_evals = 0;
    
    t0 = tspan(1); tf = tspan(end);
    t = t0;
    XA = X0;
    t_list = [t];
    X_list = [XA];
    h = h_ref;

    num_steps = 0;

    t_next = t+h;

    if t_next>tf
        t_next = tf;
        h = tf-t;
    end

    failed_steps=0;
    
    while t<tf
        

        [XB, num_evals_step, h_next, redo] = explicit_RK_variable_step(rate_func_in,t,XA,h,BT_struct,p,error_desired);
        
        if redo==false
            
            XA = XB;
            t = t_next;

            t_list(end+1) = t;
            X_list(:,end+1) = XA;

            num_steps = num_steps+1;
        else
            failed_steps=failed_steps+1;
        end


        h = h_next;
        t_next = t+h;

        if t_next>tf
            t_next = tf;
            h = tf-t;
        end

        num_evals = num_evals + num_evals_step;
        
    end

    h_avg = (tf-t0)/num_steps;

    X_list = X_list';

    step_fail_rate= failed_steps/num_steps;
end
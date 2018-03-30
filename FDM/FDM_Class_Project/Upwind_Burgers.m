function [U] = Upwind_Burgers(space_step,time_step,dx,dt_vec)
    L2_err = zeros(4,1);
    Linf_err = zeros(4,1);
    step_size = zeros(4,1);




    U = zeros(length(space_step),1);
    U_new = zeros(length(space_step),1);
    
    %set initial condition
    for i = 1:length(space_step)
        U(i) = 1/2+sin(space_step(i));
    end

    for i = 1:length(time_step)
        %Set end values not looped over
        if U(1)>0
            U_new(1) = (-dt_vec(i)/dx)*(U(1)*(U(1)-U(end)))+U(1);
        else
            U_new(1) = (-dt_vec(i)/dx)*(U(1)*(U(2)-U(1)))+U(1);
        end
        if U(end)>0
            U_new(end) = (-dt_vec(i)/dx)*(U(end)*(U(end)-U(end-1)))+U(end);
        else
            U_new(end) = (-dt_vec(i)/dx)*(U(1)*(U(1)-U(end)))+U(end);
        end
        
        %looped over values
        for j = 2:length(space_step)-1
            if U(j)>0
                U_new(j) = (-dt_vec(i)/dx)*U(j)*(U(j)-U(j-1))+U(j);
            else
                U_new(j) = (-dt_vec(i)/dx)*U(j)*(U(j+1)-U(j))+U(j);
            end
        end
        U=U_new;
    end


%     titles={'Step Number' 'L2 Error' 'L inf Error' };
%     ERR=horzcat(step_num, L2_err, Linf_err);
%     ERRMatrix=[titles; num2cell(ERR)];
%     ERRMatrix
end

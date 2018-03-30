function [U] = Lax_Friedrichs (space_step,time_step,dx,dt_vec)
    
    U = zeros(length(space_step),1);
    U_new = zeros(length(space_step),1);
    
    %set initial condition
    for i = 1:length(space_step)
        U(i) = 1/2+sin(space_step(i));
    end

    for i = 1:length(time_step)
        U_new(1) = (1/2)*(U(2)+U(end)) - (dt_vec(i)/(2*dx))*(.5*U(2)^2-.5*U(end)^2);
        U_new(end) =  (1/2)*(U(1)+U(end-1)) -(dt_vec(i)/(2*dx))*(.5*U(1)^2-.5*U(end-1)^2);
        for j = 2:length(space_step)-1
            U_new(j) = (1/2)*(U(j-1)+U(j+1)) - (dt_vec(i)/(2*dx))*(.5*U(j+1)^2-.5*U(j-1)^2);
        end
        U=U_new;
    end
end
function [U] = Godonov(space_step,time_step,dx,dt_vec)
    
    
    U = zeros(length(space_step),1);
    U_new = zeros(length(space_step),1);

            
    
    %set initial condition
    for i = 1:length(space_step)
        U(i) = 1/2+sin(space_step(i));
    end

    for i = 1:length(time_step)
        U_new(1) = U(1) - (dt_vec(i)/dx)*(flux(U(1),U(2))-flux(U(end),U(1)));
        U_new(end) = U(end) - (dt_vec(i)/dx)*(flux(U(end),U(1))-flux(U(end-1),U(end)));
        for j = 2:length(space_step)-1
            U_new(j) = U(j) - (dt_vec(i)/dx)*(flux(U(j),U(j+1))-flux(U(j-1),U(j)));
        end
        U=U_new;
    end
end
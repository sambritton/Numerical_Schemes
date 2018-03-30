function [u_newton] = Newton(f, fprime, Tol, u_0)
    
    error=1;

    tempf=f(u_0);

    while (error > Tol)  %note that dx and f need to be defined for this statement to proceed


        u_newton = u_0 - (f(u_0)/fprime(u_0));   % compute the new value of x
        error=abs(u_0-u_newton);          % compute how much x has changed since last step
        u_0 = u_newton;
        tempf=f(u_0);       % compute the new value of f(x)
        %fprintf('%3i %12.8f %12.8f %12.8f\n',count,x,dx,f)
    end
    
end
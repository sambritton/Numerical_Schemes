%FEM u_t - u_xx = f from x in [0,1] and t in [0,1]

%Set rhs function
uexact = @(x,t)((x^2-1)*(x)*(t-1)*(t));
uexact_x = @(x,t)((x^2-1)*(t-1)*t+((2*x^2)*(t-1)*t));
uexact_t = @(x,t)(x-1)*(x)*(2*t-1);
uexact_xx = @(x,t)(6*x)*t*(t-1);
f = @(x,t) ((x-1)*(x)*(2*t-1)-(6*x)*t*(t-1));


%initial condition




%Set 3 point quadrature rule
quad_num = 3;
%actual points
quad_point = [-sqrt(3/5),0,sqrt(3/5)];
%weights
quad_wgt = [5/9,8/9,5/9];


iterations = 6;
elements = zeros(iterations,1);
errorsL2 = zeros(iterations,1);
errorsH1 = zeros(iterations,1);
accuracyL2 = zeros(iterations,1);
accuracyH1 = zeros(iterations,1);


for k = 2:iterations
    %space start, end, and nodes
    a = 0;
    b = 1;
    nodes_x = (2^k)+1;
    num_elts_x = nodes_x-1;

    %time start, end, and nodes
    t_0 = 0;
    t_1 = 1;
    nodes_t = (2^(2*k))+1;
    dt = (t_1 - t_0)/(nodes_t -1);
    num_elts_t = nodes_t-1;

    space_step = linspace(a,b, nodes_x );
    time_step = linspace(t_0,t_1, nodes_t);

    %Assemble stiffness matrix, S, and Mass matrix M.
    S = zeros(nodes_x,nodes_x);
    M = zeros(nodes_x,nodes_x);
    F = zeros(nodes_x,1);
    for ex = 1:num_elts_x
        l = ex;
        r = ex + 1;
        xl = space_step(l);
        xr = space_step(r);
        dx = xr-xl;

        for q = 1:quad_num
            xq = xl + ((1+quad_point(q))/2)*dx;
            wq = quad_wgt(q)*(xr-xl)/2;

            %Basis functions
            phi_l = (xq-xl)/dx;
            phi_r = (xr-xq)/dx;
            Dphi_l = (1)/dx;
            Dphi_r = (-1)/dx;

            S(l,l) = S(l,l) + wq*(Dphi_l*Dphi_l);
            S(l,r) = S(l,r) + wq*(Dphi_l*Dphi_r);

            S(r,r) = S(r,r) + wq*(Dphi_r*Dphi_r);
            S(r,l) = S(r,l) + wq*(Dphi_r*Dphi_l);

            M(l,l) = M(l,l) + wq*(phi_l*phi_l);
            M(l,r) = M(l,r) + wq*(phi_l*phi_r);

            M(r,r) = M(r,r) + wq*(phi_r*phi_r);
            M(r,l) = M(r,l) + wq*(phi_r*phi_l);

        end
    end

    %impose boundary condition on M and S
    M(1,1:nodes_x)=0;
    M(nodes_x,1:end) = 0;
    M(2:nodes_x,1)=0;
    M(1:end,nodes_x) = 0;
    S(1,2:nodes_x)=0;
    S(nodes_x,1:end) = 0;
    S(2:nodes_x,1)=0;
    S(1:end,nodes_x) = 0;

    S(1,1) = 1;
    S(nodes_x,nodes_x) = 1;
    M(1,1) = 1;
    M(nodes_x,nodes_x) = 1;


    U=zeros(nodes_t,nodes_x);
    
    %set initial value for U using initial condition.
    for i = 1:nodes_x
        U(1,i) = sin(pi*space_step(i));
    end
    F = zeros(nodes_x,1);

    for n = 2:nodes_t

        for m = 1:nodes_x-1
            xl = space_step(m);
            xr = space_step(m+1);
            dx = xr-xl;

            for q=1:quad_num
                xq = ((1-quad_point(q))*xl + (1+quad_point(q))*xr)/2;
                wq = quad_wgt(q)*(xr-xl)/2;

                %Basis functions
                phi_l = (xq-xl)/dx;
                phi_r = (xr-xq)/dx;
                F(m) = F(m) + wq*f(xq,time_step(n))*(phi_r);

            end
        end

        %crank nicolsen
        U(n,1:end) =(M+(dt/2)*S)\((M-(dt/2)*S)*transpose(U(n-1,1:end)));
                    +transpose(F\(M+(dt/2)*S));

        % Forward Euler
        %U(n,1:end)= (M\(-dt*S*transpose(U(n-1,1:end))))+transpose(U(n-1,1:end));


        % Backward Euler
        %U(n,1:end)= (M+(dt*S))\(M*transpose(U(n-1,1:end)));


        %impose boundary conditions
        %U(n,1) = 0;
        %U(n,nodes_x) = 0;
        %reset F for next iteration in time
        F=zeros(nodes_x,1);
    end




       %Now that we have the matrix U we'll perform quadrature
    %with the psi(t)
    L2Error = 0;
    for n = 1:nodes_t-1
        t1 = time_step(n);
        t2 = time_step(n+1);
        dt = t2-t1;
        
        for m=1:nodes_x-1
            l = ex;
            r = ex + 1;
            xl = space_step(l);
            xr = space_step(r);
            dx = xr-xl;
         

            for qt=1:quad_num

                tq =((1-quad_point(qt))*t1 + (1+quad_point(qt))*t2)/2;
                wq_t = quad_wgt(qt)*(t2-t1)/2;
                
                for qx = 1:quad_num
                    xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                    wq_x = quad_wgt(qx)*(xr-xl)/2;
                    
                    wq = wq_x*dx*wq_t*dt;                    

                    phi_l = (xq-xl)/dx;
                    phi_r = (xr-xq)/dx;
                    
                    phi_b = (tq-t1)/dt;
                    phi_t = (t2-tq)/dt;
                
                    Uq = U(n,m)*phi_b + U(n+1,m)*phi_t + U(n,m)*phi_l + U(n,m+1)*phi_r;
                    eq = uexact(xq,tq);

                    L2Error = L2Error + wq*(Uq - eq)^2;
                end
            end
        end
    end
    
    H1Error = 0;
    for n=1:nodes_t-1
        t1 = time_step(n);
        t2 = time_step(n+1);
        dt = t2-t1;
        
        for m = 1:nodes_x-1
            l = ex;
            r = ex + 1;
            xl = space_step(l);
            xr = space_step(r);
            dx = xr-xl;
                
            for qt=1:quad_num
                
                tq =((1-quad_point(qt))*t1 + (1+quad_point(qt))*t2)/2;
                wq_t = quad_wgt(qt)*(t2-t1)/2;
                
                for qx=1:quad_num
                
                    xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                    wq_x = quad_wgt(qx)*(xr-xl)/2;
                    
                    %local error (don't forget to scale by dx*dt)
                    wq = wq_x*dx*wq_t*dt;
                
                    %time basis
                    Dphi_b = 1/dt;
                    Dphi_t = -1/dt;
                    
                    %space basis
                    Dphi_l = 1/dx;
                    Dphi_r = -1/dx;
                    
                    %Estimate and exact in t direction
                    Dutq = U(n,m)*Dphi_b + U(n+1,m)*Dphi_t;
                    Deq_t =  uexact_t(xq,tq);
                    
                    %Estimate and exact in x direction
                    Deq_x = uexact_x(xq,tq);
                    Duxq = U(n,m)*Dphi_l + U(n,m+1)*Dphi_r;
                    
                    
                    H1Error = H1Error + wq*((Deq_x - Duxq)^2 +(Deq_t - Dutq)^2);
                end
            end
        end
    end
    
    
    
    errorsL2(k) = sqrt(L2Error);
    errorsH1(k) = sqrt(H1Error);
    elements(k) = num_elts_x;
   
end


for i =2:iterations
    accuracyL2(i) = log2(errorsL2(i-1)/errorsL2(i));
    accuracyH1(i) = log2(errorsH1(i-1)/errorsH1(i));
end


% Display Error Matrix
titles={'Elements' 'L2 Error' 'H1 Error' 'L2 Accuracy' 'H1 Accuracy'};
ERR=horzcat(elements, errorsL2, errorsH1, accuracyL2, accuracyH1);
ERRMatrix=[titles; num2cell(ERR)];
ERRMatrix




        

    
    
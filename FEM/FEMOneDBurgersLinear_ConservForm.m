%We solve Burgers equation in conservation form. u_t*v + (u^2)_x*v = 0
%on the interval [0,1]x[0,1]
%This changes the structure of the problem. To estimate u^2 we use 
%the trial solution u^2 = sum(c_i(t)^2*phi_i(x)).
%For good introduction on the topic and numerics look at page 16 of 
%https://theses.lib.vt.edu/theses/available/etd-12052009-020403/unrestricted/LD5655.V855_1995.P844.pdf

uexact = @(x,t)(x/(1-t));
uexact_x = @(x,t) 1/(1-t);
uexact_t = @(x,t) (-x)/(1-t)^2;
f_rhs = @(x,t) 0;

%Set 3 point quadrature rule
quad_num = 3;
%actual points
quad_point = [-sqrt(3/5),0,sqrt(3/5)];
%weights
quad_wgt = [5/9,8/9,5/9];


iterations = 10;
elements = zeros(iterations,1);
errorsL2 = zeros(iterations,1);
errorsH1 = zeros(iterations,1);
accuracyL2 = zeros(iterations,1);
accuracyH1 = zeros(iterations,1);

for z=1:iterations
    %space start, end, and nodes
    a = 0;
    b = 1;
    nodes_x = (2^z)+1;
    num_elts_x = nodes_x-1;
    space_step = linspace(a,b, nodes_x );

    %time start, end, and nodes
    t_0 = 0;
    t_1 = 1/2;
    nodes_t = (2^(z))+1;
    dt = (t_1 - t_0)/(nodes_t -1);
    num_elts_t = nodes_t-1;
    time_step = linspace(t_0,t_1, nodes_t);

    %Assemble stiffness matrix, J, and Mass matrix M.
    %Notice that we will only use the stiffness matrix
    %when the viscous term is present. 
    J = zeros(nodes_x,nodes_x);
    M = zeros(nodes_x,nodes_x);
    
    for ex = 1:num_elts_x
        l = ex;
        r = ex + 1;
        xl = space_step(l);
        xr = space_step(r);
        dx = xr-xl;

        for q = 1:quad_num
            xq = xl + ((1+quad_point(q))/2)*dx;
            wq = quad_wgt(q)*(dx)/2;
            
            %Generate Basis functions
            if (2<= ex && ex <= num_elts_x)
                phi_l = (xq-xl)/dx;
                phi_r = (xr-xq)/dx;
                Dphi_l = (1)/dx;
                Dphi_r = (-1)/dx;
            elseif ex == 1
                phi_l = 0;
                phi_r = (xr-xq)/dx;
                Dphi_l = 0;
                Dphi_r = -1/dx;
            end
                
            

            %Fill matrix J and M
            J(l,l) = J(l,l) + wq*(phi_l*Dphi_l);
            J(l,r) = J(l,r) + wq*(-phi_l*Dphi_r);

            J(r,r) = J(r,r) + wq*(-phi_r*Dphi_r);
            J(r,l) = J(r,l) + wq*(phi_r*Dphi_l);

            M(l,l) = M(l,l) + wq*(phi_l*phi_l);
            M(l,r) = M(l,r) + wq*(phi_l*phi_r);

            M(r,r) = M(r,r) + wq*(phi_r*phi_r);
            M(r,l) = M(r,l) + wq*(phi_r*phi_l);

        end
    end
    
    %Make matrices nonsingular
    J(nodes_x,nodes_x) = 1/2;
    M(nodes_x,nodes_x) = dx;
    J(1,1) = 1/2;
    M(1,1) = dx;
    
    %Set Solution matrix to hold c_i values. 
    U = zeros(nodes_t,nodes_x);
    
    %Set initial condition
    for i=1:nodes_x
        U(1,i) = uexact(space_step(i),0);
    end
    
    
    %Solve
    for n = 1:num_elts_t
        t1 = time_step(n);
        t2 = time_step(n+1);
        dt = t2-t1;
            
        %Make RHS vector
        F_N_Conserv = zeros(num_elts_x+1,1);
        
        %Set values not looped over
        F_N_Conserv(1) = U(n,2)^2-U(n,1)^2;
        F_N_Conserv(end) = U(n,end)^2 - U(n,end-1)^2;
        
        %Notice that this fills in the inside of the vector only when
        %you have enough elements.
        for i = 2:num_elts_x
            F_N_Conserv(i) = U(n,i+1)^2 - U(n,i-1)^2;
        end
        
        %Solve for next row of values using B.E.
        U(n+1,1:end) = transpose(M\(-(dt/4)*F_N_Conserv)) + U(n,1:end);
        
        %boundary condition trial
        %U(n+1,1) = 0;
    end
    
    
    Q = zeros(nodes_t,nodes_x);
    for i =1:nodes_t
        for j = 1:nodes_x
            Q(i,j) = uexact(space_step(j),time_step(i));
        end
    end
    H = U-Q;

    %% L2Error
    %Now that we have the matrix U we'll perform quadrature
    %with the psi(t)
    L2Error = 0;
    t1 = time_step(nodes_t-1);
    t2 = time_step(nodes_t);
    dt = t2-t1;

    for m=1:nodes_x-1
        l = m;
        r = m + 1;
        xl = space_step(l);
        xr = space_step(r);
        dx = xr-xl;


        for qt=1:quad_num

            tq =((1-quad_point(qt))*t1 + (1+quad_point(qt))*t2)/2;
            wq_t = quad_wgt(qt)*(t2-t1)/2;

            for qx = 1:quad_num
                
                xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                wq_x = quad_wgt(qx)*(xr-xl)/2;

                %Local rectangular quadrature wgt
                wq = wq_x*wq_t;                    

                %Basis for x
                phi_l = (xq-xl)/dx;
                phi_r = (xr-xq)/dx;

                %Basis for t
                phi_b = (tq-t1)/dt;
                phi_t = (t2-tq)/dt;

                Uq = U(nodes_t-1,l)*phi_b*phi_l + U(nodes_t-1,r)*phi_b*phi_r +...
                        U(nodes_t,l)*phi_t*phi_l + U(nodes_t,r)*phi_t*phi_r;
                eq = uexact(xq,tq);

                L2Error = L2Error + wq*(Uq - eq)^2;
            end
        end
    end

%%
    H1Error = 0;
    t1 = time_step(nodes_t-1);
    t2 = time_step(nodes_t);
    dt = t2-t1;

    for m = 1:nodes_x-1
        l = m;
        r = m + 1;
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
                %How does the CFL number affect scaling? Ask!
                wq = wq_x*wq_t;

                %time basis
                phi_b = (tq-t1)/dt;
                phi_t = (t2-tq)/dt;
                Dphi_b = 1/dt;
                Dphi_t = -1/dt;

                %space basis
                phi_l = (xq-xl)/dx;
                phi_r = (xr-xq)/dx;
                Dphi_l = 1/dx;
                Dphi_r = -1/dx;

                %Estimate and exact in t direction
                Dutq = U(nodes_t-1,l)*Dphi_b*phi_l + U(nodes_t,r)*Dphi_b*phi_r+...
                        U(nodes_t,l)*Dphi_t*phi_l + U(nodes_t,r)*Dphi_t*phi_r;
                Deq_t =  uexact_t(xq,tq);

                %Estimate and exact in x direction
                Duxq = U(nodes_t-1,l)*phi_b*Dphi_l + U(nodes_t-1,r)*phi_b*Dphi_r+...
                        U(nodes_t,l)*phi_t*Dphi_l + U(nodes_t,r)*phi_t*Dphi_r ;
                Deq_x = uexact_x(xq,tq);
                


                H1Error = H1Error + wq*((Deq_x - Duxq)^2 +(Deq_t - Dutq)^2);
            end
        end
    end
    
    errorsL2(z) = sqrt(L2Error);
    errorsH1(z) = sqrt(H1Error);
    elements(z) = num_elts_x;
    
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

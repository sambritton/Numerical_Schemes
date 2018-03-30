


%Solve u_t + ku_x = 0 and u_t + ku_x = f. 
%With u(a,t) = u(b,t) = 0 and u(x,0) = initial condition.

k=-1;

%Space and Time start and end
a = 0;
b = 1;
t_0 = 0;
t_1 = 1;



uexact = @(x,t)(-1/4>t-x && t-x>-3/4);
uexact_t = @(x,t) (0);
uexact_x = @(x,t) (0);
initial_condition=@(x)(1/4<=x && x<=3/4);

f = @(x,t)(0);


%Set 3 point quadrature rule
quad_num = 3;
%actual points
quad_point = [-sqrt(3/5),0,sqrt(3/5)];
%weights
quad_wgt = [5/9,8/9,5/9];


iterations = 5;
elements = zeros(iterations,1);
errorsL2 = zeros(iterations,1);
errorsH1 = zeros(iterations,1);
accuracyL2 = zeros(iterations,1);
accuracyH1 = zeros(iterations,1);

for z=1:iterations
    %space start, end, and nodes
    nodes_x = (2^z)+1;
    num_elts_x = nodes_x-1;
    space_step = linspace(a,b, nodes_x );

    %time start, end, and nodes
    
    %question: with 2*z and 3*z there is a up/down/up in accuracy
    %with 4 it stops with the source term but not without it
    nodes_t = (2^(2*z))+1;
    dt = (t_1 - t_0)/(nodes_t -1);
    num_elts_t = nodes_t-1;
    time_step = linspace(t_0,t_1, nodes_t);

    %Assemble stiffness matrix, J, and Mass matrix M.
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
            wq = quad_wgt(q)*(xr-xl)/2;

            %Basis functions
            phi_l = (xq-xl)/dx;
            phi_r = (xr-xq)/dx;
            Dphi_l = (1)/dx;
            Dphi_r = (-1)/dx;

            %Fill matrix J and M
            J(l,l) = J(l,l) + wq*(-Dphi_l*phi_l);
            J(l,r) = J(l,r) + wq*(Dphi_l*phi_r);

            J(r,r) = J(r,r) + wq*(Dphi_r*phi_r);
            J(r,l) = J(r,l) + wq*(-Dphi_r*phi_l);

            M(l,l) = M(l,l) + wq*(phi_l*phi_l);
            M(l,r) = M(l,r) + wq*(phi_l*phi_r);

            M(r,r) = M(r,r) + wq*(phi_r*phi_r);
            M(r,l) = M(r,l) + wq*(phi_r*phi_l);

        end
    end

    %impose boundary condition on M and J
%     M(1,1:nodes_x)=0;
%     M(nodes_x,1:end) = 0;
%     M(2:nodes_x,1)=0;
%     M(1:end,nodes_x) = 0;
%     J(1,2:nodes_x)=0;
%     J(nodes_x,1:end) = 0;
%     J(2:nodes_x,1)=0;
%     J(1:end,nodes_x) = 0;
% 
%     J(1,1) = 1;
%     J(nodes_x,nodes_x) = 1;
%     M(1,1) = 1;
%     M(nodes_x,nodes_x) = 1;
% 

    U=zeros(nodes_t,nodes_x);
        
    %set initial value for U using initial condition.
    for i = 1:nodes_x
        U(1,i) = initial_condition(space_step(i));
    end
    
    %Define F for following steps
    F = zeros(nodes_x,1);

    for n = 2:nodes_t
        t1 = time_step(n-1);
        t2 = time_step(n);
        
        for m = 1:nodes_x-1
            xl = space_step(m);
            xr = space_step(m+1);
            dx = xr-xl;
            

            for qx=1:quad_num
                xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                wq_x = quad_wgt(qx)*(xr-xl)/2;
                
                for qt=1:quad_num
                    tq =((1-quad_point(qt))*t1 + (1+quad_point(qt))*t2)/2;
                    wq_t = quad_wgt(qt)*(t2-t1)/2;
                    
                    wq = wq_x*wq_t;
                    
                    %Basis functions
                    phi_b = (tq-t1)/dt;
                    phi_t = (t2-tq)/dt;
                    phi_l = (xq-xl)/dx;
                    phi_r = (xr-xq)/dx;
                    
                    
                    F(m) = F(m) + wq*f(xq,tq)*phi_r;
                    
                end
            end
        end
        %crank nicolsen
        %U(n,1:end) =(M+(dt/2)*J)\((M-(dt/2)*J)*transpose(U(n-1,1:end)));
        %           +transpose(F\(M+(dt/2)*J));
        %Backward Euler
        U(n,1:end) = (M-dt*J)\(dt*F +M*transpose(U(n-1,1:end)));
        
        %Reset F for next iteration
        F = zeros(nodes_x,1);
        
        %Impose Boundary Condition on U
        %U(n,1) = 0;
        %U(n,nodes_x) = 0;
    end       
    
    
    %Now that we have the matrix U we'll perform quadrature
    %with the psi(t)
    
    %% L2Error
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

                Uq = U(nodes_t-1,m)*phi_b + U(nodes_t,m)*phi_t +...
                        U(nodes_t,m)*phi_l + U(nodes_t,m+1)*phi_r;
                eq = uexact(xq,tq);

                L2Error = L2Error + wq*(Uq - eq)^2;
            end
        end
    end
    
    
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
                Dphi_b = 1/dt;
                Dphi_t = -1/dt;

                %space basis
                Dphi_l = 1/dx;
                Dphi_r = -1/dx;

                %Estimate and exact in t direction
                Dutq = U(nodes_t-1,m)*Dphi_b + U(nodes_t,m)*Dphi_t;
                Deq_t =  uexact_t(xq,tq);

                %Estimate and exact in x direction
                Deq_x = uexact_x(xq,tq);
                Duxq = U(nodes_t,m)*Dphi_l + U(nodes_t,m+1)*Dphi_r;


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



%Solve u_t + ku_x = 0 and u_t + ku_x = f. 
%With u(a,t) = u(b,t) = 0 and u(x,0) = uexact(x,0).

k=1;

%Zero on edges and top (1.5accuracy)
% uexact = @(x,t) x*t*(x-1)*(t-1);
% uexact_t = @(x,t) x*t*(x-1) + x*(x-1)*(t-1);
% uexact_x = @(x,t) t*(x-1)*(t-1) + x*t*(t-1);
% f = @(x,t) (x*t*(x-1) + x*(x-1)*(t-1)) + k*(t*(x-1)*(t-1) + x*t*(t-1));

%Zero on edges, not top (.5accuracy)
uexact = @(x,t) x*(x-1)*exp(t);
uexact_t = @(x,t) x*(x-1)*exp(t);
uexact_x = @(x,t) (x-1)*exp(t) + x*exp(t);
f = @(x,t) x*(x-1)*exp(t) + k*((x-1)*exp(t) + x*exp(t));

%Nonzero source term Nonzero edges boundary not set for this one
% uexact = @(x,t)(exp(-(pi)*t))*sin(pi*x);
% uexact_x = @(x,t)(exp(-(pi)*t))*pi*cos(pi*x);
% uexact_t = @(x,t)(-pi)*(exp(-(pi)*t))*sin(pi*x);
% f = @(x,t)(-pi)*(exp(-(pi)*t))*(sin(pi*x) + cos(pi*x));

% Nonzero on edges No Source Term boundary not set for this one
% uexact = @(x,t)(sin(pi*(x-k*t)));
% uexact_x = @(x,t)(pi*cos(pi*(x-k*t)));
% uexact_t = @(x,t)(-k*pi*cos(pi*(x-k*t)));
% f = @(x,t)(0);


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
    a = 0;
    b = 1;
    nodes_x = (2^z)+1;
    num_elts_x = nodes_x-1;
    space_step = linspace(a,b, nodes_x );

    %time start, end, and nodes
    t_0 = 0;
    t_1 = 1;
    nodes_t = 20*(2^(z))+1;
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
            wq = quad_wgt(q)*(dx)/2;
            
            %Basis functions
            if (2<= ex && ex <= num_elts_x-1)
                phi_l = (xq-xl)/dx;
                phi_r = (xr-xq)/dx;
                Dphi_l = (1)/dx;
                Dphi_r = (-1)/dx;
            elseif ex == 1
                phi_l = 0;
                phi_r = (xr-xq)/dx;
                Dphi_l = 0;
                Dphi_r = -1/dx;
            elseif ex == num_elts_x
                phi_l = (xq-xl)/dx;
                phi_r = 0;
                Dphi_l = 1/dx;
                Dphi_r = 0;
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
    
    %incorporate k, the constant
    %For some reason the sign of k needs to be changed 
    %if not, accuracy is completely lost.
    if k < 0
        k = -k;
    else 
        k = k;
    end
    
    J = k*J;
    
    U=zeros(nodes_t,nodes_x);
        
    %set initial value for U using initial condition.
    for i = 1:nodes_x
        U(1,i) = uexact(space_step(i),0);
    end
    
    %Define F for following steps
    F = zeros(nodes_x,1);

    for n = 1:num_elts_t
        t1 = time_step(n);
        t2 = time_step(n+1);
        dt = t2-t1;
        
        for m = 1:num_elts_x
            xl = space_step(m);
            xr = space_step(m+1);
            dx = xr-xl;
            

            for qx=1:quad_num
                xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                wq_x = quad_wgt(qx)*(dx)/2;
                
                for qt=1:quad_num
                    tq =((1-quad_point(qt))*t1 + (1+quad_point(qt))*t2)/2;
                    wq_t = quad_wgt(qt)*(dt)/2;
                    
                    wq = wq_x*wq_t;
                    
                    %Basis functions
                    if (1<= m && m <= num_elts_x)
                        phi_r = (xr-xq)/dx;
                        phi_l = (xq-xl)/dx;
                    else
                        phi_r = 0;
                        phi_l = 0;
                    end  
                    F(m) = F(m) + wq*f(xq,tq)*phi_l;
                    F(m+1) = F(m+1) + wq*f(xq,tq)*phi_r;
                    F(1) = 0;
                    F(end) = 0;
                end
            end
        end
        
        %Backward Euler
        U(n+1,1:end) = (M+dt*J)\(dt*F + M*transpose(U(n,1:end)));
        
        %Reset F for next iteration
        F = zeros(nodes_x,1);
        
        %Impose Boundary Condition on U
        %U(n+1,1) = uexact(space_step(1),time_step(n+1));
        %U(n+1,nodes_x) = uexact(space_step(end),time_step(n+1));
    end       
    
    
    %Now that we have the matrix U we'll perform quadrature
    %with the psi(t)
    
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
            wq_t = quad_wgt(qt)*(dt)/2;

            for qx = 1:quad_num
                
                xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                wq_x = quad_wgt(qx)*(dx)/2;

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
            wq_t = quad_wgt(qt)*(dt)/2;

            for qx=1:quad_num

                xq = ((1-quad_point(qx))*xl + (1+quad_point(qx))*xr)/2;
                wq_x = quad_wgt(qx)*(dx)/2;

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
                Dutq = U(nodes_t-1,l)*phi_b*Dphi_l + U(nodes_t,r)*phi_b*Dphi_r+...
                        U(nodes_t,l)*phi_t*Dphi_l + U(nodes_t,r)*phi_t*Dphi_r;
                Deq_t =  uexact_t(xq,tq);

                %Estimate and exact in x direction
                Duxq = U(nodes_t-1,l)*Dphi_b*phi_l + U(nodes_t-1,r)*Dphi_b*phi_r+...
                        U(nodes_t,l)*Dphi_t*phi_l + U(nodes_t,r)*Dphi_t*phi_r ;
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
A = zeros(nodes_t,nodes_x);
for i = 1:nodes_x
    for j = 1:nodes_t
        A(j,i) = uexact(space_step(i),time_step(j));
    end
end

% Display Error Matrix
titles={'Elements' 'L2 Error' 'H1 Error' 'L2 Accuracy' 'H1 Accuracy'};
ERR=horzcat(elements, errorsL2, errorsH1, accuracyL2, accuracyH1);
ERRMatrix=[titles; num2cell(ERR)];
ERRMatrix


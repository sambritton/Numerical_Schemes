%Solving -u_xx = f with Dirichlet boundary condition
%u(0)=u(1)=0. Quadratic elements.

f=@(x)(pi^2*sin(pi*x));
uexact=@(x)(sin(pi*x));
Duexact = @(x)(pi*cos(pi*x));


iterations = 4;

nodes = zeros(iterations,1);
elements = zeros(iterations,1);
errorsL2 = zeros(iterations,1);
errorsH1 = zeros(iterations,1);
accuracyL2 = zeros(iterations,1);
accuracyH1 = zeros(iterations,1);


%Make sure the nodes is odd and greater than three
for n=1:iterations

    num_nodes = 2*(2^n)+1;
    num_elts = (num_nodes-1)/2;
    x = linspace(0,1,num_nodes);


    %Number of quadrature points [0,1] after convertion
    quad_num = 3;
    %actual points from [-1,1]
    quad_point = [-sqrt(3/5),0,sqrt(3/5)];
    %weights
    quad_wgt = [5/9,8/9,5/9];

    %Set stiffness matrix
    A = zeros(num_nodes,num_nodes);
    g = zeros(num_nodes,1);

    %%Fill in stiffness matrix and RHS.

    for e=1:num_elts
        l = 2*e-1;
        m = 2*e;
        r = 2*e+1;
    
        xl = x(l);
        xm = x(m);
        xr = x(r);
        
        for q=1:quad_num
            xq = ((1-quad_point(q))*xl + (1+quad_point(q))*xr)/2;
            wq = quad_wgt(q)*(xr-xl)/2;
        
            %generate basis functions
            phi_l = ((xq-xm)*(xq-xr))/((xl-xm)*(xl-xr));
            phi_m = ((xq-xl)*(xq-xr))/((xm-xl)*(xm-xr));
            phi_r = ((xq-xl)*(xq-xm))/((xr-xl)*(xr-xm));
        
            Dphi_l = ((xq-xr)+(xq-xm))/((xl-xm)*(xl-xr));
            Dphi_m = ((xq-xl)+(xq-xr))/((xm-xl)*(xm-xr));
            Dphi_r = ((xq-xm)+(xq-xl))/((xr-xl)*(xr-xm));
        
            %Fill in A, stiffness matrix
            A(l,l) = A(l,l) + wq*(Dphi_l*Dphi_l);
            A(l,m) = A(l,m) + wq*(Dphi_l*Dphi_m);
            A(l,r) = A(l,r) + wq*(Dphi_l*Dphi_r);
        
            A(m,l) = A(m,l) + wq*(Dphi_m*Dphi_l);
            A(m,m) = A(m,m) + wq*(Dphi_m*Dphi_m);
            A(m,r) = A(m,r) + wq*(Dphi_m*Dphi_r);
        
            A(r,l) = A(r,l) + wq*(Dphi_r*Dphi_l);
            A(r,m) = A(r,m) + wq*(Dphi_r*Dphi_m);
            A(r,r) = A(r,r) + wq*(Dphi_r*Dphi_r);
        
            %fill g
            g(l) = g(l) + wq*phi_l*f(xq);
            g(m) = g(m) + wq*phi_m*f(xq);
            g(r) = g(r) + wq*phi_r*f(xq);
        end
    end
       

    %%Boundary conditions
    A(1,1:num_nodes) = 0;
    A(1,1) = 1;
    g(1) = 0;

    A(num_nodes,1:num_nodes) = 0;
    A(num_nodes,num_nodes) = 1;
    g(num_nodes) = 0;

    u = A\g;

    %% L2 Error
    L2Error=0;
    H1Error=0;
    for e=1:num_elts
        l = 2*e-1;
        m = 2*e;
        r = 2*e+1;
        
        xl = x(l);
        xm = x(m);
        xr = x(r);
    
        for q=1:quad_num
            xq = ((1-quad_wgt(q))*xl+(1+quad_wgt(q))*xr)/2;
            wq = quad_wgt(q)*(xr-xl)/2;
        
            %generate basis functions
            phi_l = ((xq-xm)*(xq-xr))/((xl-xm)*(xl-xr));
            phi_m = ((xq-xl)*(xq-xr))/((xm-xl)*(xm-xr));
            phi_r = ((xq-xl)*(xq-xm))/((xr-xl)*(xr-xm));
        
            Dphi_l = ((xq-xr)+(xq-xm))/((xl-xm)*(xl-xr));
            Dphi_m = ((xq-xr)+(xq-xl))/((xm-xl)*(xm-xr));
            Dphi_r = ((xq-xm)+(xq-xl))/((xr-xl)*(xr-xm));
        
            %Load approximate and true functions at the point xq
            uq = u(l)*phi_l + u(m)*phi_m + u(r)*phi_r;
            eq = uexact(xq);
        
            Duq = u(l)*Dphi_l + u(m)*Dphi_m + u(r)*Dphi_r;
            Deq = Duexact(xq);
        
            L2Error = L2Error + wq*((uq-eq)^2);
            H1Error = H1Error + wq*((Duq-Deq)^2);
        end
    end

    L2Error=sqrt(L2Error);
    H1Error=sqrt(H1Error);

    %Fill diplay vectors
    nodes(n) = num_nodes;
    elements(n) = num_elts;
    errorsL2(n) = L2Error;
    errorsH1(n) = H1Error;
end

for n=2:iterations
    accuracyL2(n) = log2((errorsL2(n-1))/(errorsL2(n)));
    accuracyH1(n) = log2(errorsH1(n-1)/errorsH1(n));
end
    
% Display Error Matrix
titles={'Elements' 'L2 Error' 'H1 Error' 'L2 Accuracy' 'H1 Accuracy'};
ERR=horzcat(elements, errorsL2, errorsH1, accuracyL2, accuracyH1);
ERRMatrix=[titles; num2cell(ERR)];
ERRMatrix
    
 
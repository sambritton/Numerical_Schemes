%Solve -u_xx = f with Dirichlet boundary conditions
%u(0)=u(1)=0. We use linear basis elements

f=@(x)(pi^2*sin(pi*x));
u=@(x)(sin(pi*x));
uPrime = @(x)(pi*cos(pi*x));
iterations = 4;

accuracy = zeros(iterations,1);
errorsl2 = zeros(iterations,1);
errorsh1 = zeros(iterations,1);
elements = zeros(iterations,1);

%Number of quadrature points [0,1] after convertion
quad_num = 3;
%actual points converted
quad_point = [(-sqrt(3/5)/2)+1/2,1/2,(sqrt(3/5)/2)+1/2];
%weights
quad_wgt = [5/18,8/18,5/18];

for q = 1:iterations
    N = 2^q;
    h = 1/N;
    x = linspace(0,1,N+1);
    

    phi = {@(x)(1-x/h);
            @(x)(x/h);};

    Dphi = {@(x)( x*0-1/h );
            @(x)( x*0+1/h ); };

    E0 = zeros(2,2);
    E2 = zeros(2,2);

    for i = 1:2
        for j = 1:2
            E2(i,j) = integral(@(X)(Dphi{j}(X).*Dphi{i}(X)),0,h);
            E0(i,j) = integral(@(X)(phi{j}(X).*phi{i}(X)),0,h);
        end;
    end;

%constructing stiffness matrix
    A = zeros(N+1,N+1);
    F = zeros(N+1,1);
    for i = 1:N
        A(i:i+1,i:i+1) = A(i:i+1,i:i+1) +E2;
        F(i:i+1) = F(i:i+1) + E0*transpose(f(x(i:i+1)));
    end;

%Impose Initial Conditions (Dirichlet) on estimate solution

    U = [0;A(2:N, 2:end-1)\F(2:end-1); 0];

    L2Error = 0;
    H1Error = 0;
    Diff = @(X)(u(X) - interp1(x,U,X));
    for i=1:quad_num
        L2Error = L2Error + quad_wgt(i)*Diff(quad_point(i));
    end;
    L2Error = sqrt(abs(L2Error));

    UPrime = zeros(N+1,1);
    for i = 1:N
        UPrime(i) = (U(i+1)-U(i))/h;
    end;

    DiffPrime = @(X)(uPrime(X) - interp1(x,UPrime,X));
    for i=1:quad_num
        H1Error = H1Error + quad_wgt(i)*DiffPrime(quad_point(i)); 
    end;
    H1Error = sqrt(abs(H1Error));

    
   
% Fill in error vectors
    
    elements(q)=N;
    errorsl2(q)=L2Error;
    errorsh1(q)=H1Error;
    
end;

for n=2:iterations
    accuracy(n)=log2(errorsl2(n-1)/errorsl2(n));
end;
% Display Error Matrix
titles={'Elements' 'L2 Error' 'H1 Error' 'Accuracy'};
ERR=horzcat(elements, errorsl2, errorsh1, accuracy);
ERRMatrix=[titles; num2cell(ERR)];
ERRMatrix
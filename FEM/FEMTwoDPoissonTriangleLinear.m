%Solve -u_xx-u_yy = f(x,y) on [0,1]x[0,1]

uexact = @(x,y) (sin(pi*x)*sin(pi*y));
uexact_x = @(x,y) (pi*cos(pi*x)*sin(pi*y));
uexact_y = @(x,y) (pi*sin(pi*x)*cos(pi*y));
f = @(x,y) ((2*pi^2)*sin(pi*x)*sin(pi*y));

iterations = 5;
order_elts = 3; %because we're using triangles????

%setting up vectors
elements = zeros(iterations,1);
errorsL2 = zeros(iterations,1);
errorsH1 = zeros(iterations,1);
accuracyL2 = zeros(iterations,1);
accuracyH1 = zeros(iterations,1);

%%
for n = 1:iterations

    %set up problem
    N = 2^n+1;
    nx = N;                 %nodes in x direction
    ny = N;                 %nodes in y direction
    
    num_nodes = nx*ny;      %Total number of nodes
    node_xy = zeros(2,num_nodes);
    
    xl = 0;
    xr = 1;
    yb = 0;
    yt = 1;
    
    %example of node coordinates with 3 rows 4 columns
    %R=3| K=9  K=10 K=11 K=12
    %R=2| K=5  K=6  K=7  K=8
    %R=1| K=1  K=2  K=3  K=4
    %   --------------------
    %     C=1  C=2  C=3  C=4
    %The following labels coordinates as pictured above
    k=0;
    for j=1:ny
        for i=1:nx
            k=k+1;
            node_xy(1,k)=((nx-i)*xl + (i-1)*xr)/(nx-1); %x coord
            node_xy(2,k)=((ny-j)*yb + (j-1)*yt)/(ny-1); %y coord
        end
    end
    
    % Element Array
    % Split eah square into two right triangles
    
    % |\  \  \ |
    % | \  \  \|
    % 5--6--7--8
    % |\  \  \ |
    % | \  \  \|
    % 1--2--3--4
    %we now break up the squares
    
    num_elts = 2*(nx-1)*(ny-1);     %total number of triangular elts
    node_elts = zeros(order_elts,num_elts);   %

    k = 0;
    for j = 1:ny-1
        for i =1:nx-1
            k=k+1;
            %bottom triangle
            node_elts(1,k) = i + (j-1)*nx;       %Bottom left
            node_elts(2,k) = (i+1) + (j-1)*nx;   %Top left
            node_elts(3,k) = i + j*nx;           %Bottom right
            
            k=k+1;
            %top triangle
            node_elts(1,k) = (i+1) + j*nx;       %Top right
            node_elts(2,k) = i + j*nx;           %Bottom right
            node_elts(3,k) = (i+1) + (j-1)*nx;   %Top left
        end
    end
    
    %% Assembly
    
    A = sparse(num_nodes,num_nodes);
    g = zeros(num_nodes,1);
    
    for e = 1:num_elts
        %we label the three vertices
        i1 = node_elts(1,e);
        i2 = node_elts(2,e);
        i3 = node_elts(3,e);
        
        %Calculate area (or delta_x)
        %               | 1 x1 y1 |
        % Area = 1/2 det| 1 x2 y2 |
        %               | 1 x3 y3 |
    
        area = (1/2)*...
            (node_xy(1,i1)*(node_xy(2,i2)-node_xy(2,i3))...
            +node_xy(1,i2)*(node_xy(2,i3)-node_xy(2,i1))...
            +node_xy(1,i3)*(node_xy(2,i1)-node_xy(2,i2)));
        
        %quadrature (1,2),(2,3),(3,1)
        for q1 = 1:3
            q2 = mod(q1,3)+1;
            
            nq1 = node_elts(q1,e);
            nq2 = node_elts(q2,e);
            
            xq = (1/2)*(node_xy(1,nq1)+node_xy(1,nq2));
            yq = (1/2)*(node_xy(2,nq1)+node_xy(2,nq2));
            wq = 1/3;
            
            %test functions
            for t_i1 = 1:order_elts
                t_i2 = mod(t_i1,3)+1;
                t_i3 = mod(t_i1+1,3)+1;
                
                nt_i1 = node_elts(t_i1,e);
                nt_i2 = node_elts(t_i2,e);
                nt_i3 = node_elts(t_i3,e);
                
                %test functions and x,y derivatives
                 phi_i = (1/2)*...
                    ((node_xy(1,nt_i3)-node_xy(1,nt_i2))*(yq-node_xy(2,nt_i2))...
                    -(node_xy(2,nt_i3)-node_xy(2,nt_i2))*(xq-node_xy(1,nt_i2)))...
                    /area;
                 phi_i_x = -(1/2)*(node_xy(2,nt_i3)-node_xy(2,nt_i2))/area;
                 phi_i_y = (1/2)*(node_xy(1,nt_i3)-node_xy(1,nt_i2))/area;
                 
                 % Fill in g
                 g(nt_i1) = g(nt_i1) + area*wq*f(xq,yq)*phi_i;
                 
                 %Basis functions
                 for t_j1 = 1:order_elts
                     t_j2 = mod(t_j1,3)+1;
                     t_j3 = mod(t_j1+1,3)+1;
                    
                     nt_j1 = node_elts(t_j1,e);
                     nt_j2 = node_elts(t_j2,e);
                     nt_j3 = node_elts(t_j3,e);
                    
                     phi_j = (1/2)*...
                         ((node_xy(1,nt_j3)-node_xy(1,nt_j2))*(yq-node_xy(2,nt_j2))...
                         -(node_xy(2,nt_j3)-node_xy(2,nt_j2))*(xq-node_xy(1,nt_j2)))...
                         /area;
                     phi_j_x = -(1/2)*(node_xy(2,nt_j3)-node_xy(2,nt_j2))/area;
                     phi_j_y =  (1/2)*(node_xy(1,nt_j3)-node_xy(1,nt_j2))/area;
                     
                     %Fill A
                     A(nt_i1,nt_j1) = A(nt_i1,nt_j1) + area*wq*(phi_i_x*phi_j_x+phi_i_y*phi_j_y);
                 end
            end
        end
    end
    
    %% Impose Boundary Conditions
    k = 0;
    for j = 1:ny
        for i = 1:nx
            k=k+1;
            if (i==1 || j==1 || i==nx || j==ny)
                A(k,1:num_nodes) = 0;
                A(k,k) = 1;
                b(k) = uexact(node_xy(1,k),node_xy(2,k));
            end
        end
    end
    u = A\g;
    u=reshape(u,nx,ny);
                
    %%Error
    L2Error = 0;
    H1Error = 0;
    
    for e = 1:num_elts
        %we label the three vertices
        i1 = node_elts(1,e);
        i2 = node_elts(2,e);
        i3 = node_elts(3,e);
        
        %Calculate area (or delta_x)
        %               | 1 x1 y1 |
        % Area = 1/2 det| 1 x2 y2 |
        %               | 1 x3 y3 |
    
        area = (1/2)*...
            (node_xy(1,i1)*(node_xy(2,i2)-node_xy(2,i3))...
            +node_xy(1,i2)*(node_xy(2,i3)-node_xy(2,i1))...
            +node_xy(1,i3)*(node_xy(2,i1)-node_xy(2,i2)));
        
        %quadrature (1,2),(2,3),(3,1)
        for q1 = 1:3
            q2 = mod(q1,3)+1;
            
            nq1 = node_elts(q1,e);
            nq2 = node_elts(q2,e);
            
            xq = (1/2)*(node_xy(1,nq1)+node_xy(1,nq2));
            yq = (1/2)*(node_xy(2,nq1)+node_xy(2,nq2));
            wq = 1/3;
        
            uh = 0;
            uh_x = 0;
            uh_y = 0;
            
            for t_j1 = 1:order_elts
                t_j2 = mod(t_j1,3)+1;
                t_j3 = mod(t_j1+1,3)+1;
                    
                nt_j1 = node_elts(t_j1,e);
                nt_j2 = node_elts(t_j2,e);
                nt_j3 = node_elts(t_j3,e);
                    
                phi_j = (1/2)*...
                     ((node_xy(1,nt_j3)-node_xy(1,nt_j2))*(yq-node_xy(2,nt_j2))...
                     -(node_xy(2,nt_j3)-node_xy(2,nt_j2))*(xq-node_xy(1,nt_j2)))...
                      /area;
                phi_j_x = -(1/2)*(node_xy(2,nt_j3)-node_xy(2,nt_j2))/area;
                phi_j_y =  (1/2)*(node_xy(1,nt_j3)-node_xy(1,nt_j2))/area;
                     
                uh = uh + u(nt_j1)*phi_j;
                uh_x = uh_x + u(nt_j1)*phi_j_x;
                uh_y = uh_y + u(nt_j1)*phi_j_y;
            end   
    
            uex = uexact(xq,yq);
            uex_x = uexact_x(xq,yq);
            uex_y = uexact_y(xq,yq);
            
            L2Error = L2Error + area*((uh-uex)^2);
            H1Error = H1Error + area*((uh_x-uex_x)^2 + (uh_y - uex_y)^2);
        end
    end
    
    L2Error = sqrt(L2Error);
    H1Error = sqrt(H1Error);
    
    elements(n) = num_elts;
    errorsL2(n) = L2Error;
    errorsH1(n) = H1Error;

end
for i=2:iterations
    accuracyL2(i) = log2(errorsL2(i-1)/errorsL2(i));
    accuracyH1(i) = log2(errorsH1(i-1)/errorsH1(i));
end
 
% Display Error Matrix
titles={'Elements' 'L2 Error' 'H1 Error' 'L2 Accuracy' 'H1 Accuracy'};
ERR=horzcat(elements, errorsL2, errorsH1, accuracyL2, accuracyH1);
ERRMatrix=[titles; num2cell(ERR)];
ERRMatrix
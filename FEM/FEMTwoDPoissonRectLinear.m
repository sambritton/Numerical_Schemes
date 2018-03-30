uexact = @(x,y)(x*(1-x)*y*(1-y));
f=@(x,y)(2*x*(x-1)+2*y*(y-1));

%To be used in H1Error estimate
uexact_x = @(x,y)((1-x)*y*(1-y)-x*y*(1-y));
uexact_y = @(x,y)(x*(1-x)*(1-y)-x*(1-x)*y);

%Size of square to be estimated
a = 0;
b = 1;

iterations = 4;

elements = zeros(iterations,1);
errorsL2 = zeros(iterations,1);
errorsH1 = zeros(iterations,1);
accuracyL2 = zeros(iterations,1);
accuracyH1 = zeros(iterations,1);

for z=1:iterations
    num_elts = 2^z;
    nodes = num_elts+1;
    side_nodes = linspace(a,b,nodes);

    allnodes = (nodes)*(nodes);

    %Number of quadrature points [0,1] after convertion
    quad_num = 3;
    %actual points converted
    quad_point = [(-sqrt(3/5)/2)+1/2,1/2,(sqrt(3/5)/2)+1/2];
    %weights
    quad_wgt = [5/18,8/18,5/18];

    %%
    %create stiffness matrix (A) and RHS g
    A = zeros(allnodes,allnodes);
    g = zeros(allnodes,1);

    %computing each element
    for ex = 1:num_elts
        %Determine width of each element
        l = ex;                 %left index
        r = ex + 1;             %right index
        xl = side_nodes(l);     %left endpoint
        xr = side_nodes(r);     %right endpoint
        dx = xr - xl;           %element width

        for ey = 1:num_elts
            %height of each square
            b = ey;                 %bottom index
            t = ey + 1;             %top index
            yb = side_nodes(b);     %bottom endpoint
            yt = side_nodes(t);     %top endpoint
            dy = yt - yb;           %element height

            % Determine node indices
            bl = (ey-1)*(nodes) + ex;
            br = (ey-1)*(nodes) + ex+1;
            tl = (ey)*(nodes) + ex;
            tr = (ey)*(nodes) + ex+1;
            %1-2,5-6,9-10,13-14,2-3,6-7
            %quadrature for each element
            %convert from [0,1]x[0,1] to [xl,xr]x[xl,xr]
            for qx = 1:3
                xq = xl + quad_point(qx)*dx;

                for qy = 1:3
                    yq = yb + quad_point(qy)*dy;
                    wq = quad_wgt(qx)*dx*quad_wgt(qy)*dy;

                    %define basis functions phi and x&y derivatives
                    phi_bl = ((xl-xq)/dx)*((yt-yq)/dy);
                    phi_bl_x = (-1/dx)*((yt-yq)/dy);
                    phi_bl_y = ((xl-xq)/dx)*(-1/dy);

                    phi_br = ((xq-xr)/dx)*((yt-yq)/dy);
                    phi_br_x = (1/dx)*((yt-yq)/dy);
                    phi_br_y = ((xq-xr)/dx)*(1/dy);

                    phi_tl = ((xl-xq)/dx)*(yq-yb)/dy;
                    phi_tl_x = (-1/dx)*((yq-yb)/dy);
                    phi_tl_y = ((xl-xq)/dx)*(-1/dy);

                    phi_tr = ((xq-xr)/dx)*((yq-yb)/dy);
                    phi_tr_x = (1/dx)*((yq-yb)/dy);
                    phi_tr_y = ((xq-xr)/dx)*(1/dy);

                    %fill in RHS values
                    g(bl) = g(bl) + wq*phi_bl*f(xq,yq);
                    g(br) = g(br) + wq*phi_br*f(xq,yq);
                    g(tl) = g(tl) + wq*phi_tl*f(xq,yq);
                    g(tr) = g(tr) + wq*phi_tr*f(xq,yq);

                    %Fill in stiffness matrix
                    A(bl,bl) = A(bl,bl) + wq*(phi_bl_x*phi_bl_x + phi_bl_y*phi_bl_y);
                    A(bl,br) = A(bl,br) + wq*(phi_bl_x*phi_br_x + phi_bl_y*phi_br_y);
                    A(bl,tl) = A(bl,tl) + wq*(phi_bl_x*phi_tl_x + phi_bl_y*phi_tl_y);
                    A(bl,tr) = A(bl,tr) + wq*(phi_bl_x*phi_tr_x + phi_bl_y*phi_tr_y);

                    A(br,bl) = A(br,bl) + wq*(phi_br_x*phi_bl_x + phi_br_y*phi_bl_y);
                    A(br,br) = A(br,br) + wq*(phi_br_x*phi_br_x + phi_br_y*phi_br_y);
                    A(br,tl) = A(br,tl) + wq*(phi_br_x*phi_tl_x + phi_br_y*phi_tl_y);
                    A(br,tr) = A(br,tr) + wq*(phi_br_x*phi_tr_x + phi_br_y*phi_tr_y);

                    A(tl,bl) = A(tl,bl) + wq*(phi_tl_x*phi_bl_x + phi_tl_y*phi_bl_y);
                    A(tl,br) = A(tl,br) + wq*(phi_tl_x*phi_br_x + phi_tl_y*phi_br_y);
                    A(tl,tl) = A(tl,tl) + wq*(phi_tl_x*phi_tl_x + phi_tl_y*phi_tl_y);
                    A(tl,tr) = A(tl,tr) + wq*(phi_tl_x*phi_tr_x + phi_tl_y*phi_tr_y);

                    A(tr,bl) = A(tr,bl) + wq*(phi_tr_x*phi_bl_x + phi_tr_y*phi_bl_y);
                    A(tr,br) = A(tr,br) + wq*(phi_tr_x*phi_br_x + phi_tr_y*phi_br_y);
                    A(tr,tl) = A(tr,tl) + wq*(phi_tr_x*phi_br_x + phi_tr_y*phi_tl_y);
                    A(tr,tr) = A(tr,tr) + wq*(phi_tr_x*phi_tr_x + phi_tr_y*phi_tr_y);

                end;
            end;
        end;
    end;

    %%
    %incorporate boundary conditions
    v = 0;
    for i=1:nodes
        for j=1:nodes
            v=v+1;
            %We need to make all edge nodes zero. If you're confused look at U 
            %After we reshape it.
            if (i==1 || j==1 || i== nodes || j == nodes)
                A(v,2:allnodes)=0;
                A(2:allnodes,v)=0;
                A(v,v)=1; %Make sure this comes after, else you'll kill it off
                g(v) = 0;
            end;

        end;

    end;

    %Solve for u
    U=A\g;  %now U is 25,1 vector.

    %we make U into a nodes by nodes matrix
    %Check that the edge of U is zeros.
    U = reshape(U,nodes,nodes);


    %%
    %solve exact solution at each 

     u_true = zeros(nodes+1,1);
     node_i = zeros(nodes+1,1);
     node_j = zeros(nodes+1,1);

     %construct node values of exact solution
     v=1;
     for i=1:nodes
         for j=1:nodes
             u_true(v)= uexact(side_nodes(i),side_nodes(j));
             node_i(v) = side_nodes(i);
             node_j(v) = side_nodes(j);
             v=v+1;
         end;
     end;


    %nodal error
    node_error=0;
    for i=1:nodes
        for j=1:nodes
            node_error=node_error+abs(uexact(i,j) - U(i,j));
        end;
    end;


    %%
    %L2Error
    L2Error=0;
    for ex = 1:num_elts
        %Determine width of each element
        l = ex;                 %left index
        r = ex + 1;             %right index
        xl = side_nodes(l);     %left endpoint
        xr = side_nodes(r);     %right endpoint
        dx = xr - xl;           %element width

        for ey = 1:num_elts
            %height of each square
            b = ey;                 %bottom index
            t = ey + 1;             %top index
            yb = side_nodes(b);     %bottom endpoint
            yt = side_nodes(t);     %top endpoint
            dy = yt - yb;           %element height

            % Determine node indices
            bl = (ey-1)*(nodes) + ex;
            br = (ey-1)*(nodes) + ex+1;
            tl = (ey)*(nodes) + ex;
            tr = (ey)*(nodes) + ex+1;

            %1-2,5-6,9-10,13-14,2-3,6-7
            %quadrature for each element
            %convert from [0,1]x[0,1] to [xl,xr]x[xl,xr]
            for qx = 1:3
                xq = xl + quad_point(qx)*dx;

                for qy = 1:3
                    yq = yb + quad_point(qy)*dy;
                    wq = quad_wgt(qx)*dx*quad_wgt(qy)*dy;

                    %define basis functions phi and x&y derivatives
                    phi_bl = ((xl-xq)/dx)*((yt-yq)/dy);
                    phi_br = ((xq-xr)/dx)*((yt-yq)/dy);
                    phi_tl = ((xl-xq)/dx)*(yq-yb)/dy;
                    phi_tr = ((xq-xr)/dx)*((yq-yb)/dy);

                    %
                    uq = U(ex,ey)*phi_bl+U(ex+1,ey)*phi_br +U(ex,ey+1)*phi_tl+U(ex+1,ey+1)*phi_tr;
                    eq = uexact(ex,ey);
                    L2Error = L2Error + wq*(uq-eq)^2;
                end;
            end;
        end;
    end;


    %%
    %H1Error
    H1Error=0;
    for ex = 1:num_elts
        %Determine width of each element
        l = ex;                 %left index
        r = ex + 1;             %right index
        xl = side_nodes(l);     %left endpoint
        xr = side_nodes(r);     %right endpoint
        dx = xr - xl;           %element width

        for ey = 1:num_elts
            %height of each square
            b = ey;                 %bottom index
            t = ey + 1;             %top index
            yb = side_nodes(b);     %bottom endpoint
            yt = side_nodes(t);     %top endpoint
            dy = yt - yb;           %element height

            % Determine node indices
            bl = (ey-1)*(nodes) + ex;
            br = (ey-1)*(nodes) + ex+1;
            tl = (ey)*(nodes) + ex;
            tr = (ey)*(nodes) + ex+1;
            %1-2,5-6,9-10,13-14,2-3,6-7
            %quadrature for each element
            %convert from [0,1]x[0,1] to [xl,xr]x[xl,xr]
            for qx = 1:3
                xq = xl + quad_point(qx)*dx;

                for qy = 1:3
                    yq = yb + quad_point(qy)*dy;
                    wq = quad_wgt(qx)*dx*quad_wgt(qy)*dy;

                    %define basis functions phi and x&y derivatives
                    phi_bl_x = (-1/dx)*((yt-yq)/dy);
                    phi_bl_y = ((xl-xq)/dx)*(-1/dy);

                    phi_br_x = (1/dx)*((yt-yq)/dy);
                    phi_br_y = ((xq-xr)/dx)*(1/dy);

                    phi_tl_x = (-1/dx)*((yq-yb)/dy);
                    phi_tl_y = ((xl-xq)/dx)*(-1/dy);

                    phi_tr_x = (1/dx)*((yq-yb)/dy);
                    phi_tr_y = ((xq-xr)/dx)*(1/dy);

                    %generate estimate solutions
                    uxq = U(ex,ey)*phi_bl_x+U(ex+1,ey)*phi_br_x+U(ex,ey+1)*phi_tl_x+U(ex+1,ey+1)*phi_tr_x;
                    uyq = U(ex,ey)*phi_bl_y+U(ex+1,ey)*phi_br_y +U(ex,ey+1)*phi_tl_y+U(ex+1,ey+1)*phi_tr_y;


                    exq = uexact_x(ex,ey);
                    eyq = uexact_y(ex,ey);

                    H1Error = H1Error + wq*((uxq-exq)^2+(uyq-eyq)^2);
                end;
            end;
        end;
    end;
    
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

                
                

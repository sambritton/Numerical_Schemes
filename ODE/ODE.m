f_1 = @(t,y) (-y);
f_2 = @(t,y) (-100*y);
f_3 = @(t,y) (((3/2)*t^2 + 2*t +1)/(y-1));

iterations = 6; 
%Set Perameters
a = 0;
b = 2; 

global num_scheme;

%Define Vectors
step_size = zeros(iterations,1);
error = zeros(iterations,1);
error_quot = zeros(iterations,1);
error_45 = zeros(iterations,1);
error_45_quot = zeros(iterations,1);
log_accuracy = zeros(iterations,1);

%Exact Function choice
prompt_function = 'Press 1 through 4 for the different functions';
fun_choice = input(prompt_function);

%exact functions and initial conditions
if fun_choice == 1
    y_0 = 1;
    y_exact = @(t) exp(-1*t);
    f = f_1;
elseif fun_choice == 2
    y_0 = 1;
    y_exact = @(t) exp(-100*t);
    f = f_2;
elseif fun_choice == 3
    y_0 = -1;
    y_exact = @(t)(1-sqrt(t^3 + 2*t^2 + 2*t +4));
    f = f_3;
elseif fun_choice == 4
    %Call file with 4th function and 4 methods rewritten
    run ODE_FUN_4
else
    fun_choice = input(prompt_function);
end

if fun_choice ~= 4
    prompt_scheme = 'Enter 1 through 4 to choose numerical schemes';
    num_scheme = input(prompt_scheme);


    %Numerical Schemes
    for j = 1:iterations
        %Step Size changing with iterations
        h = .04/(2^(j-1));
        step_size(j) = h;


        %Vector of spacial nodes
        space_step = a:h:b;
        num_nodes = (b-a)/h + 1;
        

        %initial condition and estimate solution vector
        y = zeros(num_nodes,1);
        y(1) = y_0;

        %Choice of numerical schemes decided by num_scheme var
        if num_scheme == 1
            %Solver RK_3
            for i = 1:num_nodes-1
                K1 = y(i);
                K2 = y(i) + h*(1/3)*f(space_step(i),K1);
                K3 = y(i) + h*(2/3)*f(space_step(i) + (1/3)*h,K2);
                y(i+1) = y(i) + h*((1/4)*f(space_step(i),K1)+(3/4)*f(space_step(i)+(2/3)*h,K3));
            end
            [t_45,y_45] = ode45(f,space_step,y_0);
            error(j) = abs(y_exact(b) - y(end));
            error_45(j) = abs(y_exact(b) - y_45(end));

        elseif num_scheme == 2
            %Solver Adam Bashford 4-Step
            %use exact four first terms
            y(2) = y_exact(space_step(2));
            y(3) = y_exact(space_step(3));
            y(4) = y_exact(space_step(4));
          
            %notice y(end) = y(num_nodes-4+4)
            for i = 1:num_nodes-4
                y(i+4) = y(i+3) + (h/24)*(55*f(space_step(i+3),y(i+3)) -...
                            59*f(space_step(i+2),y(i+2)) +...
                            37*f(space_step(i+1),y(i+1))-...
                            9*f(space_step(i),y(i)));
            end
            [t_45,y_45] = ode45(f,space_step,y_0);
            error(j) = abs(y_exact(b) - y(end));
            error_45(j) = abs(y_exact(b) - y_45(end));

        elseif num_scheme == 3
            %Solver Implicit Trapazoidal Rule
            for i = 1:num_nodes-1

                z = y(i);
                %Fixed point iteration
                for n = 1:100
                    z = y(i) + (1/2)*h*(f(space_step(i),y(i)) +...
                            f(space_step(i+1),z));
                end
                y(i+1) = z;
            end
            [t_45,y_45] = ode45(f,space_step,y_0);
            error(j) = abs(y_exact(b) - y(end));
            error_45(j) = abs(y_exact(b) - y_45(end));

        elseif num_scheme == 4
            %Solver BDF 3 step
            %Solve first 3 terms using RK3
            for i = 1:2
                K1 = y(i);
                K2 = y(i) + h*(1/3)*f(space_step(i),K1);
                K3 = y(i) + h*(2/3)*f(space_step(i) + (1/3)*h,K2);
                y(i+1) = y(i) + h*((1/4)*f(space_step(i),K1)+...
                            (3/4)*f(space_step(i)+(2/3)*h,K3));
            end
            %Implicit BDF 3
            for i = 1:num_nodes-3
                z = y(i+2);
                for n = 1:100
                    z = (1/11)*(18*y(i+2)-9*y(i+1)+2*y(i)+...
                            6*h*f(space_step(i+3),z));
                end
                y(i+3) = z;
            end
            [t_45,y_45] = ode45(f,space_step,y_0);
            error(j) = abs(y_exact(b)-y(end));
            error_45(j) = abs(y_exact(b) - y_45(end));

        else
            num_scheme = input(prompt_scheme);
        end

    end

    for i=1:iterations-1
        error_quot(i+1) = error(i)/error(i+1);
        log_accuracy(i+1) = log2(error(i)/error(i+1));
    end
    for i = 1:iterations
        error_45_quot(i) = error(i)/error_45(i);
    end

    % Display Error Matrix
    titles={'Step Size' 'Error' 'Error2/Error1' 'log2(e2/e1)' 'error/error45'};
    ERR=horzcat(step_size, error, error_quot, log_accuracy, error_45_quot);
    ERRMatrix=[titles; num2cell(ERR)];
    ERRMatrix
end


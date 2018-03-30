iterations = 5; 
%Set Perameters
a = 0;
b = 2; 

global num_scheme_4;

%Define Vectors
step_size = zeros(iterations,1);
error = zeros(iterations,1);
error_quot = zeros(iterations,1);
error_45 = zeros(iterations,1);
error_45_quot = zeros(iterations,1);
log_accuracy = zeros(iterations,1);

prompt_scheme = 'Enter 1 through 4 to choose numerical schemes';
num_scheme_4 = input(prompt_scheme);

syms x w z t;
y_exact = @(t) exp((-1/4)*(t^2));

w_prime = @(t,w,z) (z);
z_prime = @(t,w,z) (-1/2)*(t*z + w);
w_0 = 1;
z_0 = 0;

x = [w,z];
y_45 = @(t,x) [x(2); (-1/2)*(t*x(2) + x(1))];


%Numerical Schemes
for j = 1:iterations
    %Step Size changing with iterations
    h = .04/(2^(j-1));
    step_size(j) = h;
    
    
    %Vector of spacial nodes
    space_step = a:h:b;
    num_nodes = (b-a)/h+1;
    
    %initial condition and estimate solution vector
    w = zeros(num_nodes,1);
    w(1) = w_0;
    
    z = zeros(num_nodes,1);
    z(1) = z_0;
    
    %Choice of numerical schemes decided by num_scheme var
    if num_scheme_4 == 1
        %Solver RK_3
        for i = 1:num_nodes-1
            K1_w = w(i);
            K1_z = z(i);
            K2_w = w(i) + h*(1/3)*w_prime(space_step(i),K1_w,K1_z);
            K2_z = z(i) + h*(1/3)*z_prime(space_step(i),K1_w,K1_z);
            K3_w = w(i) + h*(2/3)*w_prime(space_step(i) + (1/3)*h,K2_w,K2_z);
            K3_z = z(i) + h*(2/3)*z_prime(space_step(i) + (1/3)*h,K2_w,K2_z);
            
            w(i+1) = w(i) + h*((1/4)*w_prime(space_step(i),K1_w,K1_z) + (3/4)*w_prime(space_step(i)+(2/3)*h,K3_w,K3_z));
            z(i+1) = z(i) + h*((1/4)*z_prime(space_step(i),K1_w,K1_z) + (3/4)*z_prime(space_step(i)+(2/3)*h,K3_w,K3_z));
        end
        [t_45,y_45_SOL] = ode45(y_45,space_step,[w_0,z_0]);
        error(j) = abs(y_exact(b) - w(end));
        error_45(j) = abs(y_exact(b) - y_45_SOL(end,1));
        
    elseif num_scheme_4 == 2
        %Solver Adam Bashford 4-Step
        %Solve first 4 terms using RK4
        for i = 1:3
            K1_w = w(i);
            K1_z = z(i);
            K2_w = w(i) + (h/2)*w_prime(space_step(i),K1_w,K1_z);
            K2_z = z(i) + (h/2)*z_prime(space_step(i),K1_w,K1_z);
            K3_w = w(i) + (h/2)*w_prime(space_step(i)+(1/2)*h,K2_w,K2_z);
            K3_z = z(i) + (h/2)*z_prime(space_step(i)+(1/2)*h,K2_w,K2_z);
            K4_w = w(i) + h*w_prime(space_step(i)+h,K3_w,K3_z);
            K4_z = z(i) + h*z_prime(space_step(i)+h,K3_w,K3_z);
            w(i+1) = w(i) + h*(1/6)*(w_prime(space_step(i),K1_w,K1_z)+...
                            2*w_prime(space_step(i)+(1/2)*h,K2_w,K2_z)+...
                            2*w_prime(space_step(i)+(1/2)*h,K3_w,K3_z)+...
                            w_prime(space_step(i)+h,K4_w,K4_z));
            z(i+1) = z(i) + h*(1/6)*(z_prime(space_step(i),K1_w,K1_z)+...
                            2*z_prime(space_step(i)+(1/2)*h,K2_w,K2_z)+...
                            2*z_prime(space_step(i)+(1/2)*h,K3_w,K3_z)+...
                            z_prime(space_step(i)+h,K4_w,K4_z));
                                   
        end
        %notice y(end) = y(num_nodes-4+4)
        for i = 1:num_nodes-4
            w(i+4) = w(i+3) + (h/24)*(55*w_prime(space_step(i+3),w(i+3),z(i+3)) -...
                        59*w_prime(space_step(i+2),w(i+2),z(i+2)) +...
                        37*w_prime(space_step(i+1),w(i+1),z(i+1))-...
                        9*w_prime(space_step(i),w(i),z(i)));
            z(i+4) = z(i+3) + (h/24)*(55*z_prime(space_step(i+3),w(i+3),z(i+3)) -...
                        59*z_prime(space_step(i+2),w(i+2),z(i+2)) +...
                        37*z_prime(space_step(i+1),w(i+1),z(i+1))-...
                        9*z_prime(space_step(i),w(i),z(i)));
        end
        [t_45,y_45_SOL] = ode45(y_45,space_step,[w_0,z_0]);
        error(j) = abs(y_exact(b) - w(end));
        error_45(j) = abs(y_exact(b) - y_45_SOL(end,1));
        
    elseif num_scheme_4 == 3
        %Solver Implicit Trapazoidal Rule
        for i = 1:num_nodes-1
            
            temp_w = w(i);
            temp_z = z(i);
            %Fixed point iteration
            for n = 1:100
                temp_w = w(i) + (1/2)*h*(w_prime(space_step(i),w(i),z(i)) +...
                        w_prime(space_step(i+1),temp_w,temp_z));
                temp_z = z(i) + (1/2)*h*(z_prime(space_step(i),w(i),z(i)) +...
                        z_prime(space_step(i+1),temp_w,temp_z));
            end
            w(i+1) = temp_w;
            z(i+1) = temp_z;
        end
        [t_45,y_45_SOL] = ode45(y_45,space_step,[w_0,z_0]);
        error(j) = abs(y_exact(b) - w(end));
        error_45(j) = abs(y_exact(b) - y_45_SOL(end,1));
        
    elseif num_scheme_4 == 4
        %Solver BDF 3 step
        %Solve first 3 terms using RK3
        for i = 1:2
            K1_w = w(i);
            K1_z = z(i);
            K2_w = w(i) + h*(1/3)*w_prime(space_step(i),K1_w,K1_z);
            K2_z = z(i) + h*(1/3)*z_prime(space_step(i),K1_w,K1_z);
            K3_w = w(i) + h*(2/3)*w_prime(space_step(i) + (1/3)*h,K2_w,K2_z);
            K3_z = z(i) + h*(2/3)*z_prime(space_step(i) + (1/3)*h,K2_w,K2_z);
            
            w(i+1) = w(i) + h*((1/4)*w_prime(space_step(i),K1_w,K1_z)+...
                        (3/4)*w_prime(space_step(i)+(2/3)*h,K3_w,K3_z));
            z(i+1) = z(i) + h*((1/4)*z_prime(space_step(i),K1_w,K1_z)+...
                        (3/4)*z_prime(space_step(i)+(2/3)*h,K3_w,K3_z));
        end
        %Implicit BDF 3
        for i = 1:num_nodes-3
            temp_w = w(i+2);
            temp_z = z(i+2);
            for n = 1:100
                temp_w = (1/11)*(18*w(i+2)-9*w(i+1)+2*w(i)+...
                        6*h*w_prime(space_step(i+3),temp_w,temp_z));
                temp_z = (1/11)*(18*z(i+2)-9*z(i+1)+2*z(i)+...
                        6*h*z_prime(space_step(i+3),temp_w,temp_z));
            end
            w(i+3) = temp_w;
            z(i+3) = temp_z;
        end
        [t_45,y_45_SOL] = ode45(y_45,space_step,[w_0,z_0]);
        error(j) = abs(y_exact(b)-w(end));
        error_45(j) = abs(y_exact(b) - y_45_SOL(end,1));
        
    else
        num_scheme_4 = input(prompt_scheme);
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

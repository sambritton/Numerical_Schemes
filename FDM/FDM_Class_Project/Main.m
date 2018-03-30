%This is the mainfile
syms x
iterations = 4;
 %Initialize Errors
 Error_Up_L2 = zeros(iterations,1);
 Error_Up_LI = zeros(iterations,1);
 Error_Lax_L2 = zeros(iterations,1);
 Error_Lax_LI = zeros(iterations,1);
 Error_Go_L2 = zeros(iterations,1);
 Error_Go_LI = zeros(iterations,1);
 Steps = zeros(iterations,1);
for n = 1:iterations
    
    N = 10*(2^(n-1));
    dx = 2*pi/N;
    dt = .5*dx; %CFL condition
    
    space_step = [0:dx:2*pi];
    time=2;
    time_step = [0:dt:time];
    time_step = [time_step,time];
    Steps(n) = N;
    dt_vec = zeros(length(time_step),1);
    for i = 1:length(time_step)-1
        dt_vec(i) = dt;
    end
    dt_vec(end) = time-time_step(end-1);
    %call newtons method
    Tol = .0001;
    U_newton = zeros(length(space_step),1);
    u_0 = .1;
    for i = 1:length(space_step)
        x=space_step(i);
        f=@(u)(1/2-u+sin(x-time*u)); 
        fprime =@(u)(-1-cos(x-time*u)*time);


        U_newton(i) = Newton(f,fprime,Tol,u_0);
        u_0 = U_newton(i);
    end
    
    %call upwind burgers
    U_up = Upwind_Burgers(space_step,time_step,dx,dt_vec);
    %step_num(n) = 10*2^(n-1);
    U_lax = Lax_Friedrichs(space_step,time_step,dx,dt_vec);
    U_go = Godonov(space_step,time_step,dx,dt_vec);
    
    %plots
    subplot(2,2,n);
    plot(space_step,U_newton,'--b',space_step,U_up, '-.r', space_step,U_go, '*', space_step,U_lax,'r')
    title(['N=' num2str((10*2^(n-1)))]);
    axis([0 2*pi -2 2])
    legend('U Newton', 'U Upwind', 'U Godonov', 'U Lax')
    
    %error estimates. 
     Error_Up_L2(n) = sqrt(sum(dx*abs(U_newton-U_up).^2));
     Error_Up_LI(n) = max(abs(U_newton-U_up));
     Error_Lax_L2(n) = sqrt(sum(dx*abs(U_newton-U_lax).^2));
     Error_Lax_LI(n) = max(abs(U_newton-U_lax));
     Error_Go_L2(n) = sqrt(sum(dx*abs(U_newton-U_go).^2));
     Error_Go_LI(n) = max(abs(U_newton-U_go));

end

Acc_Up_L2 = zeros(iterations,1);
Acc_Up_LI = zeros(iterations,1);
Ord_Up_L2 = zeros(iterations,1);
Ord_Up_LI = zeros(iterations,1);
Acc_Lax_L2 = zeros(iterations,1);
Acc_Lax_LI = zeros(iterations,1);
Ord_Lax_L2 = zeros(iterations,1);
Ord_Lax_LI = zeros(iterations,1);
Acc_Go_L2 = zeros(iterations,1);
Acc_Go_LI = zeros(iterations,1);
Ord_Go_L2 = zeros(iterations,1);
Ord_Go_LI = zeros(iterations,1);

for i = 1:iterations-1
    Acc_Up_L2(i+1) = Error_Up_L2(i)/Error_Up_L2(i+1);
    Acc_Up_LI(i+1) = Error_Up_LI(i)/Error_Up_LI(i+1);
    Ord_Up_L2(i+1) = log2(Acc_Up_L2(i+1));
    Ord_Up_LI(i+1) = log2(Acc_Up_LI(i+1));
    Acc_Lax_L2(i+1) = Error_Lax_L2(i)/Error_Lax_L2(i+1);
    Acc_Lax_LI(i+1) = Error_Lax_L2(i)/Error_Lax_LI(i+1);
    Ord_Lax_L2(i+1) = log2(Acc_Lax_L2(i+1));
    Ord_Lax_LI(i+1) = log2(Acc_Lax_LI(i+1));
    Acc_Go_L2(i+1) = Error_Go_L2(i)/Error_Go_L2(i+1);
    Acc_Go_LI(i+1) = Error_Go_LI(i)/Error_Go_LI(i+1);
    Ord_Go_L2(i+1) = log2(Acc_Go_L2(i+1));
    Ord_Go_LI(i+1) = log2(Acc_Go_LI(i+1));
end

%Tables
titles_Up={'Upwind' 'L2' 'L2 Quotient' 'L2 Order Upwind' 'LI' 'LI Quotient' 'LI Order'};
ERR_Up=horzcat(Steps, Error_Up_L2, Acc_Up_L2, Ord_Up_L2, Error_Up_LI, Acc_Up_LI, Ord_Up_LI); 
ERRMatrix_Up=[titles_Up; num2cell(ERR_Up)];
ERRMatrix_Up

titles_Lax={'Lax' 'L2' 'L2 Quotient' 'L2 Order Lax' 'LI' 'LI Quotient' 'LI Order'};
ERR_Lax=horzcat(Steps, Error_Lax_L2, Acc_Lax_L2, Ord_Lax_L2, Error_Lax_LI, Acc_Lax_LI, Ord_Lax_LI); 
ERRMatrix_Lax=[titles_Lax; num2cell(ERR_Lax)];
ERRMatrix_Lax

titles_Go={'Godenov' 'L2' 'L2 Quotient' 'L2 Order' 'LI' 'LI Quotient' 'LI Order'};
ERR_Go=horzcat(Steps, Error_Go_L2, Acc_Go_L2, Ord_Go_L2, Error_Go_LI, Acc_Go_LI, Ord_Go_LI); 
ERRMatrix_Go=[titles_Go; num2cell(ERR_Go)];
ERRMatrix_Go



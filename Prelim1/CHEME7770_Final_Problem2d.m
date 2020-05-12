%%For Problem 2d, first go to acdc.m and set a value of "s"; then, return
%%here and run code. Repeat for all three values of "s" and uncomment the
%%bottom lines for final plot.

%Initialization
x0=0;
y0=0;
z0=0;
tspan=0:0.1:50;
c0=[x0;y0;z0];

%Running ODE solver and plotting
[t,c] = ode45('acdc',tspan,c0);
x = transpose(c(:,1));
semilogy(tspan,x)
hold on

% xlabel('Time')
% ylabel('X (log scale)')
% legend('S=0.02','S=10','S=10^{5}')


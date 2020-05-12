%%For Problem 2e, first uncomment a value of "s" on this script; then, go
%%to acdc.m and uncomment "s=100." Return here and run script. Repeat for
%%each value of "s" on this script.

% s=1;
s=1000;

%Initialization
ax=3.9e-2;
ay=4.3e-3;
bx=6.1;
by=5.7;
dy=1.05;
dz=1.04;
zx=1.3e-5;
yz=11e-3;
xz=12e-2;
xy=7.9e-4;
nzx=2.32;
nxy=2;
nxz=2;
nyz=2;

%Solving for Steady State values of X,Y,Z for specified value of S
syms x y z
eqn = [(ax+bx*s)/(1+s+(z/zx)^nzx)-x == 0;...
       (ay+by*s)/(1+s+(x/xy)^nxy)-dy*y == 0;...
       1/(1+(x/xz)^nxz+(y/yz)^nyz)-dz*z == 0];
S=vpasolve(eqn,[x y z]);

x0=double(S.x);
y0=double(S.y);
z0=double(S.z);
tspan=0:0.1:100;
figure(1),clf

%Setting initial Steady State values
c1=[x0;y0;z0];
c2=1.25.*c1;
c3=0.75*c1;

%Running ODE solver for the three Steady States
[t,c] = ode45('acdc',tspan,c1);
z1 = transpose(c(:,3));
figure(1)
semilogy(tspan,z1)
hold on

[t,c] = ode45('acdc',tspan,c2);
z2 = transpose(c(:,3));
figure(1)
semilogy(tspan,z2)
hold on

[t,c] = ode45('acdc',tspan,c3);
z3 = transpose(c(:,3));
figure(1)
semilogy(tspan,z3)
hold on

%Plotting data
xlabel('Time')
ylabel('Z (log scale)')
legend('Steady state','1.25*Steady state','0.75*Steady state')

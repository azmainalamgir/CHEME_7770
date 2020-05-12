function cdot = acdc(t,c)

%%For Problem 2d, select an "s" to uncomment, then go back and run Problem
%%2d script:

% s=0.02;
% s=10;
% s=10e5;


%%For Problem 2e, uncomment the following "s", then run Problem 2e script:

% s=100;

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

%Defining variables
x=c(1);
y=c(2);
z=c(3);

%ODEs to be solved
dxdt=(ax+bx*s)/(1+s+(z/zx)^nzx)-x;
dydt=(ay+by*s)/(1+s+(x/xy)^nxy)-dy*y;
dzdt=1/(1+(x/xz)^nxz+(y/yz)^nyz)-dz*z;

cdot = [dxdt;dydt;dzdt];

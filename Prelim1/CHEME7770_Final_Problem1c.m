%Initialization
I=[0 5e-4 0.005 0.012 0.053 0.216 1.0];
I_model=0:0.01:1;
N=[19 21 41 67 86 93 93];

%Parameter values found by fitting I vs. mRNA using MATLAB Curve Toolbar
kx=0.549;
w1=0.2646;
w2=1344;
k=0.7994;
n=1.531;

%Calculating corrected mRNA levels and fitted model
mRNA=N./(6.02*10e23).*10e9.*10e8./(280*10e-7);
mRNA_model=kx.*(w1+w2.*((I_model.^n)./(k^n+I_model.^n)))./(1+w1+w2.*((I_model.^n)./(k^n+I_model.^n)));

%Plotting Data
scatter(I,mRNA)
xlabel('Inducer (mM)')
ylabel('mRNA (nmol/gDW)')
hold on
plot(I_model,mRNA_model)
legend('Data','Model')
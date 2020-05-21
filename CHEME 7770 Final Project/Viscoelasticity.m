eta = linspace(0,1,20);

v1 = zeros(1,length(eta));
v2 = zeros(1,length(eta));
v3 = zeros(1,length(eta));

for i = 1:length(eta)
    count=i
    v1(i) = MC(0.5,eta(i));
    v2(i) = MC(1,eta(i));
    v3(i) = MC(5,eta(i));
end

t = 18;
plot(eta(2:t-1),v1(2:t-1),eta(2:t-1),v2(2:t-1),eta(2:t-1),v3(2:t-1),'LineWidth',2)
ylim([-80,-20])
legend('\kappa_{sub} = 0.5 pN/nm','\kappa_{sub} = 1 pN/nm','\kappa_{sub} = 5 pN/nm')
xlabel('\eta (pN s/nm)')
ylabel('\nu_{filament} (nm/s)')
k_sub = [0.1:0.1:1,2:10,20:100];

% Parameters (obtained from Table S1)
F_m = -150;      
v_u = -120;     
n_c = 75;        
k_on = 1;        
k_off0 = 0.1;    
F_b = -2;        
k_clutch = 5;       

% Time initialization
dt = 0.01;               
tf = 100;
N = tf/dt;                 
t = 0:dt:dt*(N-1);

mean_F=zeros(size(k_sub));
mean_v=zeros(size(k_sub));
F_tot=zeros(1,N);

% Running Monte Carlo Simulation
for i = 1:length(k_sub)
    
    % Clutch initialization
    x_clutch = zeros(1,n_c);
    F_clutch = zeros(1,n_c);
    clutchstatus = zeros(1,n_c);
    k_off = zeros(1,n_c);
    
    % Substrate initialization
    x_sub = zeros(1,N);
    
    % Actin initialization
    v = zeros(1,N);
    
    % Number of attachments initialization
    n_e = sum(clutchstatus);
    n_d = n_c - n_e;
    
    for j = 1:N
        % Compute koff for all clutches
        for k = 1:n_c
            F_clutch(k) = k_clutch*(x_clutch(k)-x_sub(j));
            k_off(k) = k_off0*exp(F_clutch(k)/F_b);
        end
        
        % Calculate Rtot
        Rtot = n_d*k_on + sum(k_off.*clutchstatus);
        
        % Determine attachments and deattachments
        for k = 1:n_c
            % Clutch is attached
            if clutchstatus(k) == 1
                lim = k_off(k)/Rtot;
                check = rand(1);
                % Clutch deattaches
                if check < lim
                    clutchstatus(k) = 0;
                    x_clutch(k) = 0;
                end
                % Clutch is deattached
            elseif clutchstatus(k) == 0
                lim = k_on/Rtot;
                check = rand(1);
                % Clutch attaches
                if check < lim
                    clutchstatus(k) = 1;
                end
            end
        end
        
        % Determine clutch counts
        n_e = sum(clutchstatus);
        n_d = n_c - n_e;
        
        % Change substrate position
        x_sub(j) = (k_clutch*sum(x_clutch.*clutchstatus))/(k_sub(i) + k_clutch*n_e);
        
        % Compute actin velocity
        v(j) = v_u*(1-(k_sub(i)*x_sub(j))/F_m);
        
        % Update attached clutch position
        for k = 1:n_c
            % Clutch is attached
            if clutchstatus(k) == 1
                x_clutch(k) = x_clutch(k) + v(j)*dt;
            end
        end
        
        % Recalculate tension force for each clutch
        for k = 1:n_c
            F_clutch(k) = k_clutch*(x_clutch(k)-x_sub(j));
        end
        
        % Summing individual tension forces for all cluthes
        F_tot(j)=sum(F_clutch);
        
    end
    
    % Computing average retrograde flow and traction force
    mean_v(i)=mean(v(N-2000:N));
    mean_F(i)=mean(F_tot(N-2000:N));
    
end

yyaxis left
semilogx(k_sub,mean_v,'LineWidth',2)
yyaxis right
semilogx(k_sub,mean_F,'LineWidth',2)

xlabel('\kappa_{sub} (pN/nm)')
yyaxis left
ylabel('Mean Retrograde Flow (nm/s)')
yyaxis right
ylabel('Mean Traction Force (pN)')

%% Compliant Substrate

% Parameters (obtained from Table S1)
F_m = -150;      
v_u = -120;     
n_c = 75;        
k_on = 1;        
k_off0 = 0.1;    
F_b = -2;        
k_clutch = 5;   
k_sub = 0.01;      

% Time initialization
dt = 0.20;               
tf = 100;
N = tf/dt;                 
t = 0:dt:dt*(N-1);

% Clutch initialization
x_clutch = zeros(1,n_c);      
F_clutch = zeros(1,n_c);      
clutchstatus = zeros(1,n_c); 
k_off = zeros(1,n_c);         

% Substrate initialization
x_sub = zeros(1,N);          
x_dot = 0;                   

% Actin initialization
v = zeros(1,N);            

% Number of attachments initialization
n_e = sum(clutchstatus);     
n_d = n_c - n_e;               

figure(1), clf
figure(2), clf

% Running Monte Carlo Simulation
for i = 1:N
    % Compute koff for all clutches
    for j = 1:n_c
        F_clutch(j) = k_clutch*(x_clutch(j)-x_sub(i));
        k_off(j) = k_off0*exp(F_clutch(j)/F_b);
    end
    
    % Compute Rtot
    Rtot = n_d*k_on + sum(k_off.*clutchstatus);
    
    % Determine attachments and deattachments
    for j = 1:n_c
       % Clutch is attached 
       if clutchstatus(j) == 1          
           lim = k_off(j)/Rtot;
           check = rand(1);
           % Clutch deattaches
           if check < lim                
               clutchstatus(j) = 0;
               xclutch(j) = 0;
           end
       % Clutch is deattached    
       elseif clutchstatus(j) == 0       
           lim = k_on/Rtot;
           check = rand(1);
           % Clutch attaches
           if check < lim               
               clutchstatus(j) = 1;
           end
       end
    end
    
    % Determine clutch counts
    n_e = sum(clutchstatus);             
    n_d = n_c - n_e;                       
    
    % Change substrate position
    xold = x_sub(i);
    x_sub(i) = (k_clutch*sum(x_clutch.*clutchstatus))/(k_sub + k_clutch*n_e);
    
    % Update xdot
    x_dot = (x_sub(i)-xold)/dt;
    
    % Compute actin velocity
    v(i) = v_u*(1-(k_sub*x_sub(i))/F_m);
    
    % Update attached clutch position
    for j = 1:n_c
        % Clutch is attached
        if clutchstatus(j) == 1         
            x_clutch(j) = x_clutch(j) + v(i)*dt;
        end
    end
end

% Plotting
figure(1)
plot(t,x_sub)
xlabel('Time (s)')
ylabel('Substrate position x_{sub} (nm)')
hold on

figure(2)
plot(t,v)
xlabel('Time (s)')
ylabel('Retrograde Flow Rate \nu_{filament} (nm/s)')
hold on

%% Stiff substrate 
k_sub = 100;

% Time initialization
dt = 0.01;               
tf = 100;
N = tf/dt;                 
t = 0:dt:dt*(N-1);

x_clutch = zeros(1,n_c);      
F_clutch = zeros(1,n_c);      
clutchstatus = zeros(1,n_c);
k_off = zeros(1,n_c);         % Strain dependent detachment rate

x_sub = zeros(1,N);          % Location of substrate
x_dot = 0;                   % Strain rate of substrate

v = zeros(1,N);            % Actin speed

n_e = sum(clutchstatus);     % Number of attached clutches
n_d = n_c - n_e;               % Number of detached clutches

% Running Monte Carlo Simulation
for i = 1:N
    % Compute koff for all clutches
    for j = 1:n_c
        F_clutch(j) = k_clutch*(x_clutch(j)-x_sub(i));
        k_off(j) = k_off0*exp(F_clutch(j)/F_b);
    end
    
    % Calculate Rtot
    Rtot = n_d*k_on + sum(k_off.*clutchstatus);
    
    % Determine attachments and deattachments
    for j = 1:n_c
       % Clutch is attached 
       if clutchstatus(j) == 1          
           lim = k_off(j)/Rtot;
           check = rand(1);
           % Clutch deattaches
           if check < lim                
               clutchstatus(j) = 0;
               xclutch(j) = 0;
           end
       % Clutch is deattached    
       elseif clutchstatus(j) == 0       
           lim = k_on/Rtot;
           check = rand(1);
           % Clutch attaches
           if check < lim               
               clutchstatus(j) = 1;
           end
       end
    end
    
    % Determine clutch counts
    n_e = sum(clutchstatus);             
    n_d = n_c - n_e;                       
    
    % Change substrate position
    xold = x_sub(i);
    x_sub(i) = (k_clutch*sum(x_clutch.*clutchstatus))/(k_sub + k_clutch*n_e);
    
    % Update xdot
    x_dot = (x_sub(i)-xold)/dt;
    
    % Compute actin velocity
    v(i) = v_u*(1-(k_sub*x_sub(i))/F_m);
    
    % Update attached clutch position
    for j = 1:n_c
        % Clutch is attached
        if clutchstatus(j) == 1         
            x_clutch(j) = x_clutch(j) + v(i)*dt;
        end
    end
end

figure(1)
plot(t,x_sub)
legend ('\kappa_{sub} = 0.01 pN/nm','\kappa_{sub} = 100 pN/nm')

figure(2)
plot(t,v)
legend ('\kappa_{sub} = 0.01 pN/nm','\kappa_{sub} = 100 pN/nm')

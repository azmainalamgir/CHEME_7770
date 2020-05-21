function out =  MC(k_sub,eta)

% Parameters (obtained from Table S1)
F_m = -150;      % Total myosin pulling force [pN]
v_u = -120;      % Unloaded motor velocity [nm/s]
n_c = 75;        % Total Number of clutches []
k_on = 1;        % Clutch attachment rate constant [s^-1]
k_off0 = 0.1;    % Clutch detachment rate constant [s^-1]
F_b = -2;        % Charactheric bond rupture force [pN]
k_clutch = 5;     % Clutch spring constant [pN/nm]

% Time initialization
dt = 0.01;       % Timestep of 10 ms              
N = 100000;      % Number of timesteps           

% Clutch initialization
x_clutch = zeros(1,n_c);        % Location of clutches
F_clutch = zeros(1,n_c);        % Forces of clutches
clutchstatus = zeros(1,n_c);    % Status of clutch attachment
k_off = zeros(1,n_c);           % Detachment rate constant

% Substrate initializtion
x_sub = zeros(1,N);             % Location of substrate
x_dot = 0;                      % Substrate strain rate

% Actin initialization
v = zeros(1,N);                 % Actin speed

% Number of attachments initialization
n_e = sum(clutchstatus);        % Number of engaged clutches
n_d = n_c - n_e;                % Number of detached clutches


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
%                 x_clutch(j) = 0;
            end
        % Clutch is dettached
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
    x_old = x_sub(i);
    x_sub(i) = (k_clutch*sum(x_clutch.*clutchstatus))/(k_sub + eta*x_dot + k_clutch*n_e);
    
    % Update xdot
    x_dot = (x_sub(i)-x_old)/dt;
    
    % Compute actin velocity
    v(i) = v_u*(1-(k_sub*x_sub(i)+eta*x_dot)/F_m);
    
    % Update attached clutch location
    for j = 1:n_c
        % Clutch is attached
        if clutchstatus(j) == 1         
            x_clutch(j) = x_clutch(j) + v(i)*dt;
        end
    end
end

out = mean(v(N-20000:N));



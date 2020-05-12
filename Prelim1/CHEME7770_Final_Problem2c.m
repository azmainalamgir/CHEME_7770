%Initialization
s=0:0.01:5;
X=0:0.01:5;
Z=0:0.01:5;
S_index=[];
X_index=[];

ax=1.5;
bx=5;
zx=0.4;
nzx=2.7;
xz=1.5;
nxz=2.7;
x0=[0,0];

for ii=1:length(S)
    for jj=1:length(X)
        for kk=1:length(Z)
            sol1 = (ax+bx*S(ii))/(1+S(ii)+(Z(kk)/zx)^nzx)-X(jj);
            sol2 = 1/(1+(X(jj)/xz)^nxz)-Z(kk);
            %Determine if the two equations are approximately close to zero
            if abs(sol1) < 1e-2 && abs(sol2) < 1e-2
                %Calculate Jacobian to find stability
                J=[-1 -nzx/(zx^nzx)*Z(kk)^(nzx-1)*(ax+bx*S(ii))/(1+S(ii)+(Z(kk)/zx)^nzx)^2;...
                    -nzx/(xz^nzx)*X(jj)^(nzx-1)/(1+(X(jj)/xz)^nzx)^2 -1];
                lambda=eig(J);
                lambda_real=real(J);
                if lambda_real(1) < 0 && lambda_real(2) < 0
                        S_index = [S_index S(ii)];
                        X_index = [X_index X(jj)];
                end
            end
        end
    end
end

%Plotting data
scatter(S_index,X_index)
xlabel('S')
ylabel('X')

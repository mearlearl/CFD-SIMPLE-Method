%%% FDM HW3 %%%
clear; clc;
%%%  Set up Mixed Meshes
N = 8;
dx = 1/(N+1);
x = [dx/2:dx:1-dx/2];
xP = [0:dx:1];

%%% Boundary Conditions and Initial Guesses
Pin = 20;
Pout = 0;
c = 50 - 30*x;

% Guesses
P = linspace(Pin,Pout,N+2);
Pstar = P;
Pprime = zeros(1,N+2);
u = ones(1,N+1);
ustar = u;
uprime = zeros(1,N+1);
alphaP = 0.8;
alphau = 1-alphaP;

res = 1;
iterations = 0;

ulast = u;
Plast = P;


while res > 0.5e-8
    
    %%% Picard Linearization and forward difference in momentum equation.
    Plast = P;
    ulast = u;
    
    for i = 1:N+1

        u(i) = (Pstar(i)-Pstar(i+1))/(c(i)*dx*ustar(i));

    end
    
    %%% This new value of u solves momentum but not continuity
    ustar = u;
 
    %%% Solve for pressure Corrections using continuity
    
    A = ones(1,N+2);
    B = A;
    C = B;
    D = C;

    for i = 1:N+1
        d(i) = 1/(c(i)*dx*ustar(i));
    end    

    for i = 2:N+1
        A(i) = -d(i-1);
        B(i) = (d(i)+d(i-1));
        C(i) = -d(i);
        D(i) = (ustar(i-1)-ustar(i));
    end   

    A(1) = 0;
    A(end) = 0;
    C(1) = 0;
    C(end) = 0;
    D(1) = 0;
    D(end) = 0;

    Pprime = TDMA(A,B,C,D);

    %%% Correcting Pressure
    P = Pstar + Pprime;


    %%% Solve for corresponding velocity corrections
    for i = 1:N+1
        uprime(i)= d(i)*(Pprime(i)-Pprime(i+1));
    end

    u = ustar + uprime;

    %%% Start the cycle over again.
    res = norm(u - ulast);

    Pstar = P + alphaP*(Plast-P);%(P+Plast)/2 ;%+ (1-omega)*(P-Pstar);
    ustar = u + alphau*(ulast-u);%(u+ulast)/2 ;%+ (omega)*(u - ustar);

    iterations = iterations+1;
end

Q = sqrt(4/7);
plot(xP,-Q^2*(50*xP - (30/2)*xP.^2) + Pin,'-*',xP,P,x,Q*ones(1,N+1),'-*',x,u)
xlabel('Distance x')
ylabel('Pressure, Velocity')



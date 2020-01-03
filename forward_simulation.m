clear; close all; clc;
%% forward simulation
%% Load models and drive cycle
load Profile_adapted_2
load SS_LPV
load curve_t2
load variations
C0 = C0/10;

%%
addpath('C:\Users\s136161\Documents\ACC_2019_Dual_results')
load results_10_cell_regular_smooth

clear y x

opts             = optimset('fminsearch');
opts.TolX        = 1.e-12;

%

clear y x
start = 485;
k0 = 500;
ending = ending-(k0-start+1);
P_saved = P;
P_test = [-100*ones(1,k0) P_saved(1:ending)];
P = P_test(start:end);

% Scale drive cycle
N=10;
P = P*N;
%%
% Initial condition
x = [1;0];
x = repmat(x,N,1);
% x = repmat(x,1,20);
As = [];
Bs = [];


maxIter = 200;
u_opt = u_opt_saved(:,N*maxIter-(N-1):N*maxIter);

for i=1:length(P)
    for n=1:N
        w_opt(i,n) = u_opt(3*i-1,n);
        ubal(i,n) = u_opt(3*i,n);
    end
end

% System definition
for k = 1:length(P)
    for n = 1:N
        Af = A(x(2*n-1,k))^AB_var(n);
        Bf = B(x(2*n-1,k))*AB_var(n);
        As(2*n-1:2*n,2*n-1:2*n) = [1 0; 0 Af];
        Bs(2*n-1:2*n,1:n+1) = [1/(C0*C0_var(n)) zeros(1,n-1) 1/(C0*C0_var(n));Bf zeros(1,n-1) Bf];
        R0(n) = D(x(2*n-1,k))*D_var(n);
        Vemf(n) = EMF(x(2*n-1,k));
    end
    
    if k<length(P)
        for n=1:N
            C_left(n,k) = (x(2*n-1,k)-0.1)*C0*C0_var(n);
        end
        a = sum(R0);
        b = sum(Vemf)+repmat([0 1],1,N)*x(:,k)+R0*ubal(k,:)';
        c = -P(k);
        w(k) = (-b + sqrt(b^2-4*a*c))/(2*a);

        U = [w(k); ubal(k,:)'];
    else 
        a = sum(R0);
        b = sum(Vemf)+repmat([0 1],1,N)*x(:,k);
        c = -P(k);
        w(k) = (-b + sqrt(b^2-4*a*c))/(2*a);

        U = [w(k); zeros(N,1)];
    end
        
    x(:,k+1) = As*x(:,k) + Bs*U;
    for n=1:N
        y(n,k) = x(2*n,k) + [R0(n) zeros(1,n-1) R0(n) zeros(1,N-n)]*(U) + Vemf(n);
    end
    if max(y(:,k)<2.6)
        ending = k-1;
        break
    else
        ending = k;
    end
end

figure;hold on;grid on;
for i=1:10
    plot(y(i,:))
end

figure;hold on;
for i=1:10
    plot(x(2*i-1,:))
end

figure;hold on;
for n=1:N
   plot(C_left(n,:)) 
end

figure;hold on;
for n=1:N
   plot(ubal(:,n)) 
end


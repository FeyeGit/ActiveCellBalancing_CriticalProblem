clear; close all; clc;
%% forward simulation
%% Load models and drive cycle
% load Profile_adapted_2
load SS_LPV
load RC_1st_ARX_smooth_4
load curve_t2_smooth
load variations

figure;
plot([0:0.001:1],A([0:0.001:1]),'r')
figure;
plot([0:0.001:1],B([0:0.001:1]),'r')
figure;
plot([0:0.001:1],D([0:0.001:1]),'r')
figure;
plot([0:0.001:1],EMF([0:0.001:1]),'r')

%% Scale capacity
C0 = C0/10;

%% make smooth
A.Method = 'spline';
B.Method = 'spline';
D.Method = 'spline';

%%
addpath('C:\Users\s136161\Documents\ACC_2019_Dual_results\trust_region_in_constraint')
load Succesful_result_Edited_LPV_longer_3

clear y x

opts             = optimset('fminsearch');
opts.TolX        = 1.e-12;

%

clear y x

%%
% Initial condition
x = [1;0];
x = repmat(x,N,1);
% x = repmat(x,1,20);
As = [];
Bs = [];


maxIter = 176;
u_opt = u_opt_saved(:,N*maxIter-(N-1):N*maxIter);

for i=1:length(P)
    for n=1:N
        w_opt(i,n) = u_opt(2*i-1,n);
        ubal(i,n) = u_opt(2*i,n);
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
            C_left(n,k) = (x(2*n-1,k)-0.11)*C0*C0_var(n);
        end
        C_mean(1,k) = mean(C_left(:,k));
        for n=1:N
            C_relative(n,k) = mean(C_left(:,k))-C_left(n,k);
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
   plot(C_left(n,:).^2) 
end
plot(C_mean,'k')

figure;hold on;
for n=1:N
   plot(ubal(:,n)) 
end

figure;hold on;
for n=1:N
   plot(C_relative(n,:).^2) 
end

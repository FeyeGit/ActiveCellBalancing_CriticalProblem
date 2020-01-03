function [phi,gamma,Cm,Dm,Xiu] = MPC_predictMatrix_constant_LTV_N_lin_full(a,ending,n,C0_varied,AB_var,D_var,A,B,D,x,alpha1,alpha0)
%[Phi,Gamma] = predictMatrix(A,B,N)
%predictMatrix provides the MPC prediction matrices Phi and Gamma, given
%the discrete time system matrices A, B and the number of steps N.
[mA,~] = size(eye(2));

phi = zeros(ending*mA,mA);
gamma = zeros(ending*mA,ending*2);
Cm = zeros(ending,2*ending);
Dm = zeros(ending,2*ending);
gamma(1:2,1:2*ending) = zeros(2,2*ending);

    for ip = 1:ending
        Af = A(x(2*n-1,ip))^AB_var(n);
        Bf = B(x(2*n-1,ip))*AB_var(n);
        Ad = [1 0; 0 Af];
        Bd = [1/(C0_varied(n)) 1/(C0_varied(n));Bf Bf];
        R0(n) = D(x(2*n-1,ip))*D_var(n);
        
        if ip == 1
            phi(((ip-1)*mA+1):mA*ip,:) = eye(2);
            Xiu(ip) = 0;
        elseif ip > 1
            Af = A(x(2*n-1,ip-1))^AB_var(n);
            Ad = [1 0; alpha1(n,ip-1) Af];
            phi(((ip-1)*mA+1):mA*ip,:) = Ad*phi((((ip-1)-1)*mA+1):mA*(ip-1),:);
            Xiu(ip) = Af*(Xiu(ip-1))+alpha0(n,ip-1);
        end
        
        Cm(ip,(2*ip-1):(2*ip)) = [a(ip),1];
        Dm(ip,(2*ip-1):(2*ip)) = [R0(n),R0(n)];
        if ip == 2 
            Af = A(x(2*n-1,ip-1))^AB_var(n);
            Bf = B(x(2*n-1,ip-1))*AB_var(n);
            Bd = [1/(C0_varied(n)) 1/(C0_varied(n));Bf Bf];
            gamma(2*ip-1:2*ip,1:2*(ip-1)) = Bd;
        elseif ip > 2
            Af = A(x(2*n-1,ip-1))^AB_var(n);
            Bf = B(x(2*n-1,ip-1))*AB_var(n);
            Ad = [1 0; alpha1(n,ip-1) Af];
            Bd = [1/(C0_varied(n)) 1/(C0_varied(n));Bf Bf];
            gamma(2*ip-1:2*ip,1:2*(ip-1)) = [Ad*gamma(2*(ip-1)-1:2*(ip-1),1:2*(ip-2)) Bd];   
        end
    end
end
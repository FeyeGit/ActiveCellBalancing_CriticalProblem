function [omega,psi,xi] = predictMatrix_variable_full_lin(nu,mu,gamma_lambda,Rb,lambda,w_hat,N)
%[Phi,Gamma] = predictMatrix(A,B,N)
%predictMatrix provides the MPC prediction matrices Phi and Gamma, given
%the discrete time system matrices A, B and the number of steps N.
omega = zeros(N,2*N);
psi = zeros(2*N,2*N);
xi = zeros(1,2*N);
    for ip = 1:N
        omega(ip,(2*ip-1):(2*ip)) = [nu(ip),mu(ip)];
        psi((2*ip-1):(2*ip),(2*ip-1):(2*ip)) = [gamma_lambda/2,0;0,(mu(ip))*Rb];
        xi(1,(2*ip-1):(2*ip)) = [lambda(ip)-gamma_lambda*w_hat(ip),0];
    end
end
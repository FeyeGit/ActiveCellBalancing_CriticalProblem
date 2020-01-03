function [out] = MPC_acceptable_bal_currents(x_opt,N,u_opt)
    out=[u_opt-x_opt(2:N+1)]'*[u_opt-x_opt(2:N+1)];
%     out=sum(abs([u_opt-x_opt(2:N+1)]));
end


function [out] = I_complete_par_simplified_forward(x_opt,C0_varied,x,N)
    for n=1:N
        x_next(n) = x(2*n-1)+1/C0_varied(n)*(x_opt(1)+x_opt(n+1));
        C_left_next(n) = (x_next(n)-0.1)*C0_varied(n);
    end
    out = mean(abs(C_left_next-mean(C_left_next)));
    if out < 1
        out = 0;
    end
end


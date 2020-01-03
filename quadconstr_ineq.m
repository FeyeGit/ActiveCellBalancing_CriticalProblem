function [y,yeq] = quadconstr_ineq(x,H_cons_1,k_cons_1,d_cons_1,H_cons_2,k_cons_2,d_cons_2)
y = [];
yeq = zeros(1,1);
yeq(1) = 1/2*x'*H_cons_1*x + k_cons_1'*x + d_cons_1;
y(1) = 1/2*x'*H_cons_2*x + k_cons_2'*x + d_cons_2;

% grady = [];
% if nargout > 2
%     gradyeq = zeros(length(x),2);
%     gradyeq(:,1) = H_cons_1*x + k_cons_1;
%     gradyeq(:,2) = H_cons_2*x + k_cons_2;
% end
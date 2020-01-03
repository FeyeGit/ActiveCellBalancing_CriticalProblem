clear;close all;clc;
%% Simulatie results from matlab online
% addpath('D:\Matlab online')
% addpath('D:\Documents\Graduation\Results\Dual')
% load('D:\Matlab online\result_delta_bound.mat')
% % load('D:\Documents\Graduation\Results\Dual\results_bal')
% load('D:\Documents\Graduation\Results\Dual\results_10_cell_gamma_s.mat')
% load('results_adapted_cost_hybrid')
% load('results_LTV_full_lin_1')

addpath('C:\Users\s136161\Documents\ACC_2019_Dual_results\trust_region_in_constraint')
load edited_LPV_linear_piecewise_start_P_600

maxIter=300;
% N =2;
green = [12 195 82] ./ 255;
darkblue = [1 17 181] ./ 255;
red = [255 0 0]./255;
color = [darkblue;red;green];

teal = [18 150 155] ./ 255;
lightgreen = [94 250 81] ./ 255;
green = [12 195 82] ./ 255;
lightblue = [8 180 238] ./ 255;
darkblue = [1 17 181] ./ 255;
yellow = [251 250 48] ./ 255;
peach = [251 111 66] ./ 255;
brown = [224 211 162]./255;
color = [teal;lightgreen;green;lightblue;darkblue;yellow;peach;[29 11 5]./255;[249 37 138]./255;[231 8 1]./255];

figure;hold on
for n=1:N
    plot(y(n,1:ending),'Color',color(n,:))
end
for n=1:N
    plot(y_saved(n,1:ending),'Color',color(n,:),'LineStyle','--')
end

figure;hold on;
for n=1:N
    plot(x(2*n-1,1:end-1),'Color',color(n,:))
end
for n=1:N
    plot(x_saved(2*n-1,1:end-1),'Color',color(n,:),'LineStyle','--')
end
% for n=1:N
%     plot(x_base(2*n-1,1:end),'Color',color(n,:),'LineStyle','--')
% end

% figure;hold on;
% for n=1:N
%     plot(x_saved(2*n-1,1:end-1),'Color',color(n,:))
% end
% for n=1:N
%     plot(x_base(2*n-1,1:end),'Color',color(n,:),'LineStyle','--')
% end

% figure;hold on
% for n=1:N
%     plot(y_saved(n,1:end),'Color',color(n,:))
% end
% for n=1:N
%     plot(y_base(n,1:end),'Color',color(n,:),'LineStyle','--')
% end
% 
% figure;hold on; grid;
% for n=1:N
%     plot(y_matrix_saved(N*(maxIter-1)+n,:),'Color',color(n,:),'LineStyle','--')
% end

% figure;hold on;
% for n=1:N
%     plot(x(2*n-1,1:end-1),'Color',color(n,:))
% end

figure;hold on;
for n=1:N
    plot(x_saved(2*N*(maxIter-1)+2*n-1,1:end-1))
end

u_opt = u_opt_saved(:,N*maxIter-(N-1):N*maxIter);

for i=1:ending
    for n=1:N
        w_opt(i,n) = u_opt(3*i-1,n);
        ubal(i,n) = u_opt(3*i,n);
    end
end

figure;hold on;
for n=1:N
    plot(w_opt(:,n),'Color',color(n,:))
end
for n=1:N
    plot(w(1:end),'k','LineStyle','--')
end

figure;hold on;
for n=1:N
    plot(ubal(:,n),'Color',color(n,:))
end
% 
% Rb = 0.5;
% u1_opt = ubal(:,1);
% u2_opt = ubal(:,2);
% 
% for k=1:ending
%     efficiency(k) = (y(2,k)*abs(u2_opt(k)))/(y(1,k)*abs(u1_opt(k)));
% end
% 
% for k=1:ending
%     efficiency_2(k) = (y(1,k)*abs(u1_opt(k))-Rb*(u1_opt(k)^2+u2_opt(k)^2))/(y(1,k)*abs(u1_opt(k)));
% end

% figure
% plot(efficiency)


%% Multiplier evolution
mult = 254/(N*maxIter);
figure
hold on
subplot(3,1,1);hold on;grid
for iter = 1:maxIter
    color1 = [255-mult*iter 0 0] ./ 255;
    for n=1:N
        plot(lambda_saved(N*iter-(N-n),:),'Color',color1)
    end
%     pause
end
subplot(3,1,2);hold on;grid
for iter = 1:N*maxIter
    color1 = [mult*(N*maxIter)+1-mult*iter 17 181] ./ 255;
    plot(mu_saved(iter,:),'Color',color1)
%     pause
end
subplot(3,1,3);hold on;grid
for iter = 1:N*maxIter
    color1 = [mult*(N*maxIter)+1-mult*iter 17 181] ./ 255;
    plot(nu_saved(iter,:),'Color',color1)
%     pause
end

figure
hold on
subplot(3,1,1);hold on;grid
for iter = maxIter
    for n=1:N
        plot(lambda_saved(N*iter-(N-n),:),'Color',color(n,:))
    end
end
subplot(3,1,2);hold on;grid
for iter = N*maxIter
    color1 = [mult*(N*maxIter)+1-mult*iter 17 181] ./ 255;
    plot(mu_saved(iter,:),'Color',color1)
end
subplot(3,1,3);hold on;grid
for iter = N*maxIter
    color1 = [mult*(N*maxIter)+1-mult*iter 17 181] ./ 255;
    plot(nu_saved(iter,:),'Color',color1)
end

%% Error evolution
figure
hold on
subplot(3,1,1);hold on;grid
for iter = 1:maxIter
    color1 = [255-mult*iter 0 0] ./ 255;
    for n=1:N
        plot(error_lambda_saved(N*iter-(N-n),:),'Color',color1)
    end
%     pause
end
subplot(3,1,2);hold on;grid
for iter = 1:N*maxIter
    color1 = [mult*N*maxIter+1-mult*iter 17 181] ./ 255;
    plot(mu_error_saved(iter,:),'Color',color1)
%     pause
end
subplot(3,1,3);hold on;grid
for iter = 1:N*maxIter
    color1 = [mult*N*maxIter+1-mult*iter 17 181] ./ 255;
    plot(nu_error_saved(iter,:),'Color',color1)
%     pause
end

figure
hold on
subplot(3,1,1);hold on;grid
for iter = maxIter
    for n=1:N
        plot(error_lambda_saved(N*iter-(N-n),:),'Color',color(n,:))
    end
end
subplot(3,1,2);hold on;grid
for iter = N*maxIter
    color1 = [mult*N*maxIter+1-mult*iter 17 181] ./ 255;
    plot(mu_error_saved(iter,:),'Color',color1)
end
subplot(3,1,3);hold on;grid
for iter = 2994
    color1 = [mult*N*maxIter+1-mult*iter 17 181] ./ 255;
    plot(nu_error_saved(iter,:),'Color',color1)
%     pause
end

%% Check for valid solution
succesfull = [];
for iter=1:N*maxIter
    if sum(mu_error_saved(iter,:)<=0)==numel(mu_error_saved(iter,:)) && sum(nu_error_saved(iter,:)<=0)==numel(nu_error_saved(iter,:))
        succesfull = [succesfull iter];
    end
    mu_succes(iter) = sum(mu_error_saved(iter,:)<=0);
    nu_succes(iter) = sum(nu_error_saved(iter,:)<=0);
end

figure;hold on;grid
plot(mu_succes)
plot(nu_succes)

%% summed errors
for iter=1:maxIter
    for n=1:N
        summed_lambda_error(iter,n) = sum(error_lambda_saved(iter+(n-1),:));
    end
end
for iter=1:N*maxIter
    summed_mu_error(iter) = sum(mu_error_saved(iter,:));
    summed_nu_error(iter) = sum(nu_error_saved(iter,:));
end

figure
hold on
subplot(3,2,1);hold on;grid
for n=1:N
    plot(summed_lambda_error(:,n),'Color',color(n,:))
end
subplot(3,2,3);hold on;grid
plot(summed_mu_error,'Color',color(1,:))
subplot(3,2,5);hold on;grid
plot(summed_nu_error,'Color',color(1,:))
subplot(3,2,2);hold on;grid;
for n=1:N
    plot([maxIter-49:maxIter],summed_lambda_error(maxIter-49:maxIter,n),'Color',color(n,:))
end
xlim([maxIter-25 maxIter]);
subplot(3,2,4);hold on;grid;
xlim([N*maxIter-25 N*maxIter]);
plot([N*maxIter-49:N*maxIter],summed_mu_error(1,N*maxIter-49:N*maxIter),'Color',color(1,:))
subplot(3,2,6);hold on;grid;
plot([N*maxIter-49:N*maxIter],summed_nu_error(1,N*maxIter-49:N*maxIter),'Color',color(1,:))
xlim([N*maxIter-25 N*maxIter]);

for iter=1:maxIter
   total_lambda_summed(iter) = sum(summed_lambda_error(iter,:)); 
   total_mu(iter) = sum(summed_mu_error(1,N*(iter-1)+1:N*iter));
   total_nu(iter) = sum(summed_nu_error(1,N*(iter-1)+1:N*iter));
end

figure;
subplot(1,2,1)
hold on
plot(total_lambda_summed+total_mu+total_nu,'k')
plot(total_lambda_summed,'r')
plot(total_mu,'b')
plot(total_nu,'m')
legend('total','\lambda','\mu','\nu')
subplot(1,2,2)
hold on
plot([maxIter-49:maxIter],total_lambda_summed(1,maxIter-49:maxIter)+total_mu(1,maxIter-49:maxIter)+total_nu(1,maxIter-49:maxIter),'k')
plot([maxIter-49:maxIter],total_lambda_summed(1,maxIter-49:maxIter),'r')
plot([maxIter-49:maxIter],total_mu(1,maxIter-49:maxIter),'b')
plot([maxIter-49:maxIter],total_nu(1,maxIter-49:maxIter),'m')


figure;
hold on
plot(summed_mu_error+summed_nu_error)

%% multiplied with lagrangian multipliers

for iter=1:maxIter
    for n=1:N
        summed_lambda_error(iter,n) = sum(lambda_saved(iter+(n-1),:)*error_lambda_saved(iter+(n-1),:)');
    end
end
for iter=1:N*maxIter
    summed_mu_error(iter) = sum(mu_saved(iter,:)*mu_error_saved(iter,:)');
    summed_nu_error(iter) = sum(nu_saved(iter,:)*nu_error_saved(iter,:)');
end

for iter=1:maxIter
   total_lambda_summed(iter) = sum(summed_lambda_error(iter,:)); 
   total_mu(iter) = sum(summed_mu_error(1,N*(iter-1)+1:N*iter));
   total_nu(iter) = sum(summed_nu_error(1,N*(iter-1)+1:N*iter));
end

pt = 10;

addpath('C:\Users\s136161\OneDrive - TU Eindhoven\PhD\DISC\Figure requirements')
addpath('C:\Users\s136161\OneDrive - TU Eindhoven\PhD\DISC\Figure requirements\legendflex')
addpath('C:\Users\s136161\OneDrive - TU Eindhoven\PhD\DISC\Figure requirements\setgetpos_V1.2')

figure;
subplot(1,2,1)
grid
hold on
h1=plot(total_lambda_summed+total_mu+total_nu,'k');
h2=plot(total_lambda_summed,'r');
h3=plot(total_mu,'b');
h4=plot(total_nu,'m');
legend('total','\lambda','\mu','\nu')
subplot(1,2,2)
grid
hold on
plot([maxIter-49:maxIter],total_lambda_summed(1,maxIter-49:maxIter)+total_mu(1,maxIter-49:maxIter)+total_nu(1,maxIter-49:maxIter),'k')
plot([maxIter-49:maxIter],total_lambda_summed(1,maxIter-49:maxIter),'r')
plot([maxIter-49:maxIter],total_mu(1,maxIter-49:maxIter),'b')
plot([maxIter-49:maxIter],total_nu(1,maxIter-49:maxIter),'m')
legendflex([h1 h2 h3 h4],{'total','$\lambda$','$\mu$','$\nu$'},...
    'Interpreter','latex','ncol',1,'nrow',4,'anchor',[6 6],'buffer',[0 0],'padding',[0,0,3],'xscale',[0.7]);
mlf2pdf(gcf,'fig_smooth')

figure;
plot(total_lambda_summed)

figure;
hold on
plot(summed_mu_error+summed_nu_error)

%% Critical evolution
figure;plot(critics)

figure;hold on;
for iter=1:maxIter
    plot(critical_saved(:,N*iter-(N-1)))
end

%% Convergence check
% figure;hold on;grid;
% subplot(2,1,1);hold on;
% for iter = 1:3
%     for n=1:N
%         plot(x_saved(2*N*iter-(4-((2*n))),:),'Color',color(n,:))
%         plot(x_2_matrix_saved(2*(iter-1)+n,:),'Color',color(n,:),'LineStyle','--')
%     end
%     pause
% end
% % % subplot(2,1,2);hold on;grid;
% figure;hold on; grid;
% for iter = maxIter:maxIter
%     for n=1:N
%         plot(y_saved(N*iter-(1+(1-n)),:),'Color',color(n,:))
%         plot(y_matrix_saved(2*(iter-1)+n,:),'Color',color(n,:),'LineStyle','--')
%     end
% %     pause
% end

%% Figures to be exported
pt = 12;
sizeX = 12;
sizeY = 10;

% save results_20_errors_rough summed_mu_error summed_nu_error 

% save results_19_figures summed_mu_error summed_nu_error 

% save results_3_cell_ubal ubal

%{
name = 'Figures/multipliers';

mult = 254/maxIter;
multipliers = figure;
hold on
subplot(3,1,1);hold on;grid
for iter = 1:maxIter
    color1 = [mult*maxIter+1-mult*iter 17 181] ./ 255;
    plot(lambda_saved(2*iter-1,:),'Color',color1)
    color2 = [255-mult*iter 0 0] ./ 255;
    plot(lambda_saved(2*iter,:),'Color',color2)
end
ylabel('$\lambda$ [-]','FontUnits','points','interpreter','latex','FontSize',pt);
set(gca,'TickLabelInterpreter','latex','FontSize',pt)
hfig = gcf();
set(hfig,'Units','centimeters','NumberTitle','off');grid on;
subplot(3,1,2);hold on;grid
for iter = 1:maxIter
    color = [mult*maxIter+1-mult*iter 17 181] ./ 255;
    plot(mu_saved(iter,:),'Color',color)
end
ylabel('$\mu$ [-]','FontUnits','points','interpreter','latex','FontSize',pt);
set(gca,'TickLabelInterpreter','latex','FontSize',pt)
hfig = gcf();
set(hfig,'Units','centimeters','NumberTitle','off');grid on;
subplot(3,1,3);hold on;grid
for iter = 1:maxIter
    color = [mult*maxIter+1-mult*iter 17 181] ./ 255;
    plot(nu_saved(iter,:),'Color',color)
end
ylabel('$\nu$ [-]','FontUnits','points','interpreter','latex','FontSize',pt);
xlabel('Time $k$ [s]','FontUnits','points','interpreter','latex','FontSize',pt);
set(gca,'TickLabelInterpreter','latex','FontSize',pt)
hfig = gcf();
set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% name = replaceBetween('Figures/','/',fig_name);
print(hfig,name,'-painters','-depsc')

%%
name = 'Figures/errors';

mult = 254/maxIter;
multipliers = figure;
hold on
subplot(3,1,1);hold on;grid
for iter = 1:maxIter
    color1 = [mult*maxIter+1-mult*iter 17 181] ./ 255;
    plot(error_lambda_saved(2*iter-1,:),'Color',color1)
    color2 = [255-mult*iter 0 0] ./ 255;
    plot(error_lambda_saved(2*iter,:),'Color',color2)
end
ylabel('$\lambda$ error [A]','FontUnits','points','interpreter','latex','FontSize',pt);
set(gca,'TickLabelInterpreter','latex','FontSize',pt)
hfig = gcf();
set(hfig,'Units','centimeters','NumberTitle','off');grid on;
subplot(3,1,2);hold on;grid
for iter = 1:maxIter
    color = [mult*maxIter+1-mult*iter 17 181] ./ 255;
    plot(mu_error_saved(iter,:),'Color',color)
end
ylabel('$\mu$ error [W]','FontUnits','points','interpreter','latex','FontSize',pt);
set(gca,'TickLabelInterpreter','latex','FontSize',pt)
hfig = gcf();
set(hfig,'Units','centimeters','NumberTitle','off');grid on;
subplot(3,1,3);hold on;grid
for iter = 1:maxIter
    color = [mult*maxIter+1-mult*iter 17 181] ./ 255;
    plot(nu_error_saved(iter,:),'Color',color)
end
ylabel('$\nu$ error [W]','FontUnits','points','interpreter','latex','FontSize',pt);
xlabel('Time $k$ [s]','FontUnits','points','interpreter','latex','FontSize',pt);
set(gca,'TickLabelInterpreter','latex','FontSize',pt)
hfig = gcf();
set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% name = replaceBetween('Figures/','/',fig_name);
print(hfig,name,'-painters','-depsc')
%}

%%
% name = 'Figures/errors_max_LTV_6';
% 
% mult = 254/maxIter;
% multipliers = figure;
% hold on
% subplot(3,1,1);hold on;grid
% for iter = 1:maxIter
%     color1 = [mult*maxIter+1-mult*iter 17 181] ./ 255;
%     plot(error_lambda_saved(2*iter-1,:),'Color',color1)
%     color2 = [255-mult*iter 0 0] ./ 255;
%     plot(error_lambda_saved(2*iter,:),'Color',color2)
% end
% ylabel('$\lambda$ error [A]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% subplot(3,1,2);hold on;grid
% for iter = 1:maxIter
%     color1 = [mult*maxIter+1-mult*iter 17 181] ./ 255;
%     plot(mu_error_saved(iter,:),'Color',color1)
% end
% ylabel('$\mu$ error [W]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% subplot(3,1,3);hold on;grid
% for iter = 1:maxIter
%     color1 = [mult*maxIter+1-mult*iter 17 181] ./ 255;
%     plot(nu_error_saved(iter,:),'Color',color1)
% end
% ylabel('$\nu$ error [W]','FontUnits','points','interpreter','latex','FontSize',pt);
% xlabel('Time $k$ [s]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% % pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% % name = replaceBetween('Figures/','/',fig_name);
% print(hfig,name,'-painters','-depsc')
% 
% name = 'Figures/bal_current_LTV_6';
% figure;hold on
% for n=1:N
%     plot(ubal(:,n),'Color',color(n,:))
% end
% ylabel('$u$ balancing current [A]','FontUnits','points','interpreter','latex','FontSize',pt);
% xlabel('Time $k$ [s]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% % pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% % name = replaceBetween('Figures/','/',fig_name);
% print(hfig,name,'-painters','-depsc')
% 
% 
% % 
% name = 'Figures/results_summed_LTV_6';
% figure
% hold on
% subplot(3,2,1);hold on;grid
% for n=1:N
%     plot(summed_lambda_error(:,n),'Color',color(n,:))
% end
% xlim([0 maxIter]);
% ylabel('$\lambda$ error [A]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% subplot(3,2,3);hold on;grid
% plot(summed_mu_error,'Color',color(1,:))
% xlim([0 maxIter]);
% ylabel('$\mu$ error [W]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% subplot(3,2,5);hold on;grid
% plot(summed_nu_error,'Color',color(1,:))
% xlim([0 maxIter]);
% ylabel('$\nu$ error [W]','FontUnits','points','interpreter','latex','FontSize',pt);
% xlabel('Time $k$ [s]','FontUnits','points','interpreter','latex','FontSize',pt);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% xlabel('Iteration $i$ [-]','FontUnits','points','interpreter','latex','FontSize',pt);
% subplot(3,2,2);hold on;grid;
% for n=1:N
%     plot([maxIter-50:maxIter],summed_lambda_error(maxIter-50:maxIter,n),'Color',color(n,:))
% end
% xlim([maxIter-50 maxIter]);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% subplot(3,2,4);hold on;grid;
% xlim([maxIter-50 maxIter]);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% plot([maxIter-50:maxIter],summed_mu_error(1,maxIter-50:maxIter),'Color',color(1,:))
% subplot(3,2,6);hold on;grid;
% plot([maxIter-50:maxIter],summed_nu_error(1,maxIter-50:maxIter),'Color',color(1,:))
% xlim([maxIter-50 maxIter]);
% set(gca,'TickLabelInterpreter','latex','FontSize',pt)
% hfig = gcf();
% set(hfig,'Units','centimeters','NumberTitle','off');grid on;
% xlabel('Iteration $i$ [-]','FontUnits','points','interpreter','latex','FontSize',pt);
% pos = get(hfig,'position');set(hfig,'position',[pos(1:2),sizeX,sizeY]);
% % name = replaceBetween('Figures/','/',fig_name);
% print(hfig,name,'-painters','-depsc')

clear;close all;clc;
%% Model differences
load RC_1st_ARX_smooth_4
load curve_t2_smooth
s = [0:0.001:1];
A.Method = 'spline';
B.Method = 'spline';
D.Method = 'spline';

figure;
subplot(2,4,1);
plot(s,EMF(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(EMF(s)) max(EMF(s))],'r');
plot([0.3 0.3+1e-8],[min(EMF(s)) max(EMF(s))],'r')
xlabel('s')
ylabel('EMF')
subplot(2,4,2);
plot(s,A(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(A(s)) max(A(s))],'r');
plot([0.3 0.3+1e-8],[min(A(s)) max(A(s))],'r')
xlabel('s')
ylabel('A')
subplot(2,4,3);
plot(s,B(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(B(s)) max(B(s))],'r');
plot([0.3 0.3+1e-8],[min(B(s)) max(B(s))],'r')
xlabel('s')
ylabel('B')
subplot(2,4,4);
plot(s,D(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(D(s)) max(D(s))],'r');
plot([0.3 0.3+1e-8],[min(D(s)) max(D(s))],'r')
xlabel('s')
ylabel('D')
load SS_LPV
Anew = griddedInterpolant([0 0.05 0.09 A.GridVectors{1,1}],[0.75 0.78 0.8 A.Values]);
Bnew = griddedInterpolant([0 0.05 0.09 B.GridVectors{1,1}],[30e-3 20e-3 10e-3 B.Values]);
Dnew = griddedInterpolant([0 0.05 0.09 D.GridVectors{1,1}],[0.031 0.028 0.025 D.Values]);
A=Anew;
B=Bnew;
D=Dnew;
A.Method = 'spline';
B.Method = 'spline';
D.Method = 'spline';
A.ExtrapolationMethod = 'nearest';
subplot(2,4,5);
plot(s,EMF(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(EMF(s)) max(EMF(s))],'r');
plot([0.3 0.3+1e-8],[min(EMF(s)) max(EMF(s))],'r')
xlabel('s')
ylabel('EMF')
subplot(2,4,6);
plot(s,A(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(A(s)) max(A(s))],'r');
plot([0.3 0.3+1e-8],[min(A(s)) max(A(s))],'r')
xlabel('s')
ylabel('A')
subplot(2,4,7);
plot(s,B(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(B(s)) max(B(s))],'r');
plot([0.3 0.3+1e-8],[min(B(s)) max(B(s))],'r')
xlabel('s')
ylabel('B')
subplot(2,4,8);
plot(s,D(s));grid;hold on;
plot([0.1 0.1+1e-8],[min(D(s)) max(D(s))],'r');
plot([0.3 0.3+1e-8],[min(D(s)) max(D(s))],'r')
xlabel('s')
ylabel('D')

pause

load RC_1st_ARX_smooth_4
load curve_t2_smooth
s = [0:0.001:1];
A.Method = 'spline';
B.Method = 'spline';
D.Method = 'spline';

grid1 = [linspace(0.1,0.13,3) 0.15 0.18];
values = [D(grid1(1)) D(grid1(2))+0.003 D(grid1(3)) D(grid1(4)) D(grid1(5))];
GridVectors = sort([D.GridVectors{1,1};grid1']);
Values = [D.Values(1:7);values';D.Values(8:end)];
Dnew = griddedInterpolant(GridVectors,Values);
Dnew.Method = 'spline';
D=Dnew;
subplot(2,4,4)
plot(s,D(s),'m');

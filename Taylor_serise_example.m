% this program is used for showing a satalite about the earth by lagrange apprach
close all; clear all; clc;
%% constants
mu=398600;
%% Intial condition
R{1}(1,:)=[7000,0,0];
V{1}(1,:)=[0,8.2663,0];
%% solution parameter
final=900;   % final time of soution
step=6;            % time step
%% exact solution
theta=linspace(0,0.3181*pi,length(0:step:final));
r_ex=334824e4/mu./(1+0.2*cos(theta));
%% solution by first term of lagrange apprach
tl1=0;
for l1=1:length(0:step:final)
    Rl1(l1,1:3)=R{1}(1,:)+V{1}(1,:)*tl1;
    tl1=tl1+step;
end
%% solution by two terms of lagrange apprach
tl2=0;
for l2=1:length(0:step:final)
    Rl2(l2,1:3)=R{1}(1,:)+V{1}(1,:)*tl2+(-mu/(norm(R{1}(1,:)))^3*R{1}(1,:))*tl2^2/2;
    tl2=tl2+step;
end
%% solution by three terms of lagrange apprach
tl3=0;
for l3=1:length(0:step:final)
    Rl3(l3,1:3)=R{1}(1,:)+V{1}(1,:)*tl3+0.5*-mu/(norm(R{1}(1,:)))^3*R{1}(1,:)*tl3^2+(-mu*V{1}(1,:)/norm(R{1}(1,:))^3+3*mu*dot(R{1}(1,:),V{1}(1,:))/norm(R{1}(1,:))^5*R{1}(1,:))*tl3^3/6;
    tl3=tl3+step;
end
%% solution by four terms of lagrange apprach
tl4=0;
for l4=1:length(0:step:final)
    Rl4(l4,1:3)=(1-mu/2/(norm(R{1}(1,:)))^3*tl4^2+mu/2/(norm(R{1}(1,:)))^5*dot(R{1}(1,:),V{1}(1,:))*tl4^3+mu/24*(-2*mu/(norm(R{1}(1,:)))^6+3*(norm(V{1}(1,:)))^2/(norm(R{1}(1,:)))^5-15*(dot(R{1}(1,:),V{1}(1,:)))^2/(norm(R{1}(1,:)))^7)*tl4^4)*R{1}(1,:)+(tl4-mu/6/(norm(R{1}(1,:)))^3*tl4^3+mu/4*dot(R{1}(1,:),V{1}(1,:))*tl4^4/(norm(R{1}(1,:)))^5)*V{1}(1,:);
    tl4=tl4+step;
end
%% plotting
%--------------------------------------------------------------------------------------------------------------------------------------------------------
figure(1);
set(gcf,'Color','w');
hold all;
plot(0:step:final,sqrt((r_ex.*cos(theta)).^2+(r_ex.*sin(theta)).^2),'LineWidth',2)
plot(0:step:final,sqrt(Rl1(:,1).^2+Rl1(:,2).^2),0:step:final,sqrt(Rl2(:,1).^2+Rl2(:,2).^2),0:step:final,sqrt(Rl3(:,1).^2+Rl3(:,2).^2),0:step:final,sqrt(Rl4(:,1).^2+Rl4(:,2).^2),'LineWidth',2);
grid on;
legend('exact solution','solution of one term of Taylor series','solution of two terms of Taylor series','solution of three terms of Taylor series','solution of four terms of Taylor series');
xlabel('time','Fontsize',18);
ylabel('r','Fontsize',18);
%--------------------------------------------------------------------------------------------------------------------------------------------------------
%% error plotting
%--------------------------------------------------------------------------------------------------------------------------------------------------------
er1=abs(sqrt(Rl1(:,1).^2+Rl1(:,2).^2)-r_ex');
er2=abs(sqrt(Rl2(:,1).^2+Rl2(:,2).^2)-r_ex');
er3=abs(sqrt(Rl3(:,1).^2+Rl3(:,2).^2)-r_ex');
er4=abs(sqrt(Rl4(:,1).^2+Rl4(:,2).^2)-r_ex');
figure(2);
hold all;
set(gcf,'Color','w');
plot(0:step:final,er1,0:step:final,er2,0:step:final,er3,0:step:final,er4);
grid on;
legend('error of one term of Taylor series','error of two terms of Taylor series','error of three terms of Taylor series','error of four terms of Taylor series');
xlabel('time','Fontsize',18);
ylabel('error','Fontsize',18);
%--------------------------------------------------------------------------------------------------------------------------------------------------------
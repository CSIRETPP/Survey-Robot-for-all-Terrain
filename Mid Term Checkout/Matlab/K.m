clear;
clc;

%%建模思路->
%%          首先对轮腿机器人进行物理建模，因为涉及多个状态向量的处理，因此采取LQR的控制算法进行整体控制


syms theta phi x x_b N N_f T T_p M N_M P_M L_M P t
syms theta_dot x_dot phi_dot theta_ddot x_ddot phi_ddot
syms x_b_dot x_b_ddot

%%参数设定
g=9.81;%重力加速度

M=2;
L_M=0.01;
%驱动轮

mw=1;
R=0.075;
Iw=1*0.075^2;
IM=0.1;

%腿部
%大腿
l_active_leg=0.25;
m_active_leg=0.2;

%小腿
l_slave_leg=0.15;
m_slave_leg=0.2;

%摆杆重心到轮轴距离
L=0.2;

%整条腿
mp=(m_active_leg+m_slave_leg)*2+0.18;
Ip=mp*L^2;

%机体重心到转轴的距离
l=0.015;

%%经典力学分析

% x_b_dot=diff(x,t);
% x_b_ddot=diff(x_b_dot,t);
% theta_dot=diff(theta,t);
% theta_ddot=diff(theta_dot,t);

%求解常量
% equA=N-N_M-mp*diff(x+L*sin(theta),t,2)==0
% equB=P-P_M-mp*g==mp*diff(L*cos(theta),2)
% equC=N_M==M*diff((x+(L+L_M)*sin(theta)-l*sin(phi)),t,2)
% equD=P_M-M*g==M*diff((L+L_M)*cos(theta)+l*cos(phi),t,2)

equA = subs(N-N_M-mp*diff(x+L*sin(theta),t,2)==0, [diff(x,t,2), diff(theta,t,2), diff(phi,t,2)], [x_ddot, theta_ddot, phi_ddot])
equB = subs(P-P_M-mp*g==mp*diff(L*cos(theta),2), [diff(x,t,2), diff(theta,t,2), diff(phi,t,2)], [x_ddot, theta_ddot, phi_ddot])
equC = subs(N_M==M*diff((x+(L+L_M)*sin(theta)-l*sin(phi)),t,2), [diff(x,t,2), diff(theta,t,2), diff(phi,t,2)], [x_ddot, theta_ddot, phi_ddot])
equD = subs(P_M-M*g==M*diff((L+L_M)*cos(theta)+l*cos(phi),t,2), [diff(x,t,2), diff(theta,t,2), diff(phi,t,2)], [x_ddot, theta_ddot, phi_ddot])
%%消去中间变量 P，N，P_M,N_M
Equ=[equA,equB,equC,equD];
constant=solve(Equ,[P,N,P_M,N_M])


equ1=diff(x,t,2)==(T-N*R)/(Iw/R+mw*R);
equ2=Ip*diff(theta,t,2)==(P*L+P_M*L_M)*sin(theta)-(N*L+N_M*L_M)*cos(theta)-T+T_p;
equ3=IM*diff(phi,t,2)==T_p+N_M*l*cos(phi)+P_M*l*sin(phi);

equ1=subs(equ1,[P,N,P_M,N_M,diff(x,t,2),diff(phi,t,2),diff(theta,t,2)],[constant.P,constant.N,constant.P_M,constant.N_M,x_b_ddot,phi_ddot,theta_ddot])
equ2=subs(equ2,[P,N,P_M,N_M,diff(x,t,2),diff(phi,t,2),diff(theta,t,2)],[constant.P,constant.N,constant.P_M,constant.N_M,x_b_ddot,phi_ddot,theta_ddot])
equ3=subs(equ3,[P,N,P_M,N_M,diff(x,t,2),diff(phi,t,2),diff(theta,t,2)],[constant.P,constant.N,constant.P_M,constant.N_M,x_b_ddot,phi_ddot,theta_ddot])

%%求雅克比矩阵
U=[T,T_p].';
model_sol = solve([equ1,equ2,equ3],[theta_ddot,x_ddot,phi_ddot])
X = [theta,theta_dot,x,x_dot,phi,phi_dot].';
dX = [theta_dot,simplify(model_sol.theta_ddot),x_dot,simplify(model_sol.x_ddot),phi_dot,simplify(model_sol.phi_ddot)].'
    

A=subs(jacobian(dX,X),[theta theta_dot x_b_dot phi phi_dot],zeros(1,5))
B=subs(jacobian(dX,U),[theta theta_dot x_b_dot phi phi_dot],zeros(1,5))

%%LQR
%%设置权重
L=[1,1,20,300,5,10];
Q=[3,10];
K=lqr(A,B,Q,R)


L0=0.15;%腿%质心到转轴距离
%%变长
for i=1:20
    L0=L0+0.005; % 10mm线性化一次
    leglen(i)=L_var*2;
    trans_A=double(subs(A_sym,[L L_M],[L0 L0]));
    trans_B=double(subs(B_sym,[L L_M],[L0 L0]));
    KK=lqrd(trans_A,trans_B,Q_cost,R_cost,0.001);
    KK_t=KK.';
    K(i,:)=KK_t(:);
end
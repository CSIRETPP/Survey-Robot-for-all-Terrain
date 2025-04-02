clc;
clear all;
close all;

%l0变化范围
l0s=0.04:0.01:0.14;
Ks=zeros(2,6,length(l0s));%存放不同的l0对应的增益矩阵K


for step=1:length(l0s)

syms theta(t) phi(t) x(t) theta_dot phi_dot x_dot theta_dddot phi_ddot x_ddot
syms T Tp N P NM PM Nf 
syms xdd thetadd phidd

%%重力加速度
g=9.81;

%%驱动轮半径
R=0.08;

%%摆杆重心到驱动轮轴距离
L=l0s(step)/2;

%%摆杆重心到机体转轴距离
Lm=l0s(step)/2;

%%机体重心到其转轴的距离
l=0.05;

%%驱动轮转子质量
mw=0.5;

%%摆杆质量
mp=0.2;

%%机体质量
M=0.1;

%%驱动轮转子转动惯量
Iw=mw*R*R;                  %(算)

%%摆杆绕质心转动惯量
Ip=1/3*mp*((2*Lm)^2);       %（算）

%%机体绕质心转动惯量
Im=M*l*l;                   %（算）

%%求导关系
theta_dot=diff(theta,t);
theta_ddot=diff(theta,t,2);
phi_dot=diff(phi,t);
phi_ddot=diff(phi,t,2);
x_dot=diff(x,t);
x_ddot=diff(x,t,2);

%%通过对系统进行物理建模，我们可以得到中间量如下：

%求出PM
PM=M*g+M*diff((L+Lm)*cos(theta)+l*cos(phi),t,2);

%求出NM
NM=M*diff(x+(L+Lm)*sin(theta)-l*sin(phi),t,2);

%求出N
N=NM+mp*diff(x+L*sin(theta),t,2);

%求出P
P=PM+mp*g+mp*diff(L*cos(theta),t,2);



%%状态空间变量

U=[T,Tp].';


syms xdd thetadd phidd;

eq1=xdd==(T-N*R)/(Iw/R+mw*R);

eq2=Ip*thetadd==(P*L+PM*Lm)*sin(theta)-(N*L+NM*Lm)*cos(theta)-T+Tp;

eq3=Im*phidd==Tp+NM*l*cos(phi)+PM*l*sin(phi);

%解方程得到三个未知数的解theta x phi
[xdd,thetadd,phidd]=solve([eq1,eq2,eq3],[xdd,thetadd,phidd]);

x_ddot=xdd;
theta_ddot=thetadd;
phi_ddot=phidd;


%%状态空间模型
X=[theta,theta_dot,x,x_dot,phi,phi_dot].';

dX=[theta_dot,theta_ddot,x_dot,x_ddot,phi_dot,phi_ddot].';


JA = jacobian(dX,[theta,theta_dot,x,x_dot,phi,phi_dot]);
JB = jacobian(dX,[T,Tp]);
A=vpa(subs(JA,[theta,theta_dot,x,x_dot,phi,phi_dot],[0,0,x,0,0,0]));
B=vpa(subs(JB,[theta,theta_dot,x,x_dot,phi,phi_dot],[0,0,x,0,0,0]));


%离散化处理
[G,H]=c2d(eval(A),eval(B),0.005);

%%LQR
%设置权重
Q=diag([1,1,10,1000,2,3]);
R=diag([4,1]);



%计算增益反馈矩阵K
Ks(:,:,step)=dlqr(G,H,Q,R);



end


%对K的每个元素关于l0的拟合
K=sym('K',[2,6]);
syms L0;
for x=1:2
    for y=1:6
        p=polyfit(l0s,reshape(Ks(x,y,:),1,length(l0s)),3);
        K(x,y)=p(1)*L0^3+p(2)*L0^2+p(3)*L0+p(4);
    
    end
end


%matlabFunction()
matlabFunction(K,'File','lqr_K');





%%VMC求解髋关节电机力矩T1 T2
syms phi1(t) phi2(t) phi3(t) phi4(t) L5 phi0 L0 T F
syms phi1_dot phi4_dot
syms L1 L2 L3 L4

%%位置坐标表示
XB=L1*cos(phi1);
YB=L1*sin(phi1);

XD=L5+L4*cos(phi4);
YD=L4*sin(phi4);

XC=XB+L2*cos(phi2);
YC=YB+L2*sin(phi2);

XB_dot=diff(XB,t);
YB_dot=diff(YB,t);
XC_dot=diff(XC,t);
YC_dot=diff(YC,t);
XD_dot=diff(XD,t);
YD_dot=diff(YD,t);


%%速度求解
phi2_dot=(XD_dot-XB_dot)*cos(phi3)+(YD_dot-YB_dot)*sin(phi3)/(L2*sin(phi3-phi2));

x_dot_C = subs(x_dot_C,diff(phi2,t),phi_dot_2);
x_dot_C = subs(x_dot_C,[diff(phi1,t),diff(phi4,t)],[phi_dot_1,phi_dot_4]);
y_dot_C = subs(y_dot_C,diff(phi2,t),phi_dot_2);
y_dot_C = subs(y_dot_C,[diff(phi1,t),diff(phi4,t)],[phi_dot_1,phi_dot_4]);

%%运动映射+及坐标转换
X_dot_c=[x_dot_C,y_dot_C]';
q_dot_c=[phi1_dot,phi4_dot]';

J=simplify(jocobian(X_dot_C,q_dot_C));
R=[cos(phi0-pi/2),-sin(phi0-pi/2);sin(phi0-pi/2),cos(phi0-pi/2)];
Axis_Transformation=[0,-1/L0;1,0];


[T1,T2]'=J'*R*Axis_Transformation*[F,Tp]';
clear;
clc;

syms phi1(t) phi2(t) phi3(t) phi4(t) phi0 L0 T F
syms phi1_dot phi4_dot
syms L1 L2 L3 L4 L5 

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

XC_dot = subs(XC_dot,diff(phi2,t),phi2_dot);
XC_dot = subs(XC_dot,[diff(phi1,t),diff(phi4,t)],[phi1_dot,phi4_dot]);
YC_dot = subs(YC_dot,diff(phi2,t),phi2_dot);
YC_dot = subs(YC_dot,[diff(phi1,t),diff(phi4,t)],[phi1_dot,phi4_dot]);

%%运动映射+及坐标转换
x_dot=[XC_dot;YC_dot];
q_dot=[phi1_dot;phi4_dot];
x_dot = simplify(collect(x_dot,q_dot));
J=simplify(jacobian(x_dot,q_dot));
R=[cos(phi0-pi/2),-sin(phi0-pi/2);sin(phi0-pi/2),cos(phi0-pi/2)];
M=[0,-1/L0;1,0];

T = simplify(J.'*R*M)
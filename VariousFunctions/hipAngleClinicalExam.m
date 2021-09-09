function [ANG]=hipAngleClinicalExam(X,U,X2)
% [ANG]=hipAngleClinicalExam(X,U,X2)
% X..........[3x1] position of the marker heel wrt to the hip center of rotation
% U..........[3x1] axis of rotation
% X2.........[1x1] final value of the X value of the heel marker 
% define the anglo of the hip to be applied in the clinical exam tests of the knee
% it assumes that we want to have the marker on the foot at the same X value of the pelvis marker (expressed in pelvic frame -> X==0)
syms        ang
u       =   sym('u',[3,1]);
x1      =   sym('x1',[3,1]);
R       =   [cos(ang)+(u(1)^2)*(1-cos(ang))               u(1)*u(2)*(1-cos(ang))-u(3)*sin(ang)          u(1)*u(3)*(1-cos(ang))+u(2)*sin(ang)
             u(2)*u(1)*(1-cos(ang))+u(3)*sin(ang)         cos(ang)+(u(2)^2)*(1-cos(ang))                u(2)*u(3)*(1-cos(ang))-u(1)*sin(ang)
             u(3)*u(1)*(1-cos(ang))-u(2)*sin(ang)         u(3)*u(2)*(1-cos(ang))+u(1)*sin(ang)          cos(ang)+(u(3)^2)*(1-cos(ang))      ];
     
F       =   subs(solve(R(1,:)*x1==X2,ang),x1,X);
F       =   symfun(F,u);
ANG     =   F(U(1),U(2),U(3));
ANG     =   real(double(ANG));
I       =   abs(ANG)==min(abs(ANG));
ANG     =   ANG(I);

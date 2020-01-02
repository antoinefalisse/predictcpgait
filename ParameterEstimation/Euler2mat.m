%% [mat]=Euler2mat(r)
% mat:  rotation matrix - Rx*Ry*Rz
% r:    [rx, ry, rz] radiansts

function [OUT]=Euler2mat(r,ordine)
D=      size(r);
if (D(2)==3)&&(D(1)~=3)
    r=      r';
end
OUT=    zeros(3,3,size(r,2));
for q=1:size(r,2)
    X=      r(1,q);
    Y=      r(2,q);
    Z=      r(3,q);
    
    cx=     cos(X);
    sx=     sin(X);
    cy=     cos(Y);
    sy=     sin(Y);
    cz=     cos(Z);
    sz=     sin(Z);
    
    %     Rx=[1 0 0;0 cx -sx;0 sx cx];
    %     Ry=[cy 0 sy;0 1 0;-sy 0 cy];
    %     Rz=[cz -sz 0;sz cz 0;0 0 1];
    
    if nargin==1
        %         OUT(:,:,q)=Rx*Ry*Rz;
        OUT(:,:,q)=[cy*cz,                       -cy*sz,                 sy
                    cx*sz + cz*sx*sy,             cx*cz - sx*sy*sz,     -cy*sx
                    sx*sz - cx*cz*sy,             cz*sx + cx*sy*sz,      cx*cy];
        
    else
        if nargin==2
            Rx=[1 0 0;      0 cx -sx;   0 sx cx];
            Ry=[cy 0 sy;    0 1 0;      -sy 0 cy];
            Rz=[cz -sz 0;   sz cz 0;    0 0 1];
            switch ordine
                case 'zxy'
                    OUT(:,:,q)=     (Rz*Rx)*Ry;
                case 'ZXY'
                    OUT(:,:,q)=     Ry*(Rx*Rz);
                case 'zyx'
                    OUT(:,:,q)=     (Rz*Ry)*Rx;
                case 'ZYX'
                    OUT(:,:,q)=     Rx*(Ry*Rz);
                case 'xyz'
                    OUT(:,:,q)=     (Rx*Ry)*Rz;
                case 'XYZ'
                    OUT(:,:,q)=     Rz*(Ry*Rx);
            end
        end
    end
end
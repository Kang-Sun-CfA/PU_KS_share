function SG = F_2D_SG_rotate(xmesh,ymesh,x0,y0,FWHMx,FWHMy,k,Angle)
% two dimensional rotated super gaussian, only one shape parameter

% written by Kang Sun on 2017/12/09


wx = FWHMx/2/(log(2))^(1/k);
wy = FWHMy/2/(log(2))^(1/k);

if length(Angle) > 1
    xym1 = [xmesh(:)'-x0; ymesh(:)'-y0];
    xym2 = Angle*xym1;
    xx = xym2(1,:);
    yy = xym2(2,:);
    SG = exp(-abs(sqrt(xx.^2+(wx/wy*yy).^2)/wx).^k);
    SG = reshape(SG(:),size(xmesh,1),size(xmesh,2));
else
    rotation_matrix = [cos(Angle), -sin(Angle);
        sin(Angle),  cos(Angle)];
    xym1 = [xmesh(:)'-x0; ymesh(:)'-y0];
    xym2 = rotation_matrix*xym1;
    xx = xym2(1,:);
    yy = xym2(2,:);
    SG = exp(-abs(sqrt(xx.^2+(wx/wy*yy).^2)/wx).^k);
    SG = reshape(SG(:),size(xmesh,1),size(xmesh,2));
end
function SG = F_2D_SG_rotate_KY(xmeshr,ymeshr,xmesh,ymesh,x0,y0,FWHMx,FWHMy,k,Angle)
% super gaussian for circular fov. use Kai Yang's discretization.
% updated from F_2D_SG_rotate.m by Kang Sun on 2018/07/17

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
    SGc = exp(-abs(sqrt(xx.^2+(wx/wy*yy).^2)/wx).^k);
    SGc = reshape(SGc(:),size(xmesh,1),size(xmesh,2));
    
    xym1 = [xmeshr(:)'-x0; ymeshr(:)'-y0];
    xym2 = rotation_matrix*xym1;
    xx = xym2(1,:);
    yy = xym2(2,:);
    SGr = exp(-abs(sqrt(xx.^2+(wx/wy*yy).^2)/wx).^k);
    SGr = reshape(SGr(:),size(xmeshr,1),size(xmeshr,2));
    SG = (2*SGc+conv2(SGr,ones(2,2),'valid'))/6;
end
% make sure all mass is equal
meshsize = size(xmesh);
if ismember(0,meshsize) || (meshsize(1) == 1 && meshsize(2) == 1)
    return
elseif ismember(1,meshsize)
    SG = SG/(sum(SG(:))*(xmesh(2)-xmesh(1))*(ymesh(2)-ymesh(1)))...
    *FWHMx/2*FWHMy/2*pi;
else
SG = SG/(sum(SG(:))*(xmesh(1,2)-xmesh(1,1))*(ymesh(2,1)-ymesh(1,1)))...
    *FWHMx/2*FWHMy/2*pi;
end
function SG = F_fit_IASI_SRF(coeff,inp)

% forward function to fit the tabulated IASI SRF by a rotated super
% Gaussian

% Written by Kang Sun on 2017/12/08

ymesh = inp.ymesh;
zmesh = inp.zmesh;
F_SG = inp.F_SG;

w = coeff(1);
k = coeff(2);

yyzz = coeff(3);
scale = coeff(4);
SG = F_SG(sqrt((ymesh/yyzz).^2+(zmesh).^2),w,k);
SG = scale*SG(:);
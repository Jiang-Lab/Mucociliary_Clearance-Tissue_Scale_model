function Lc_r = length_rod(phi,L)
Lc_r = -0.027732*phi.^3 + 0.445795*phi.^2 -2.285569*phi + 4.240372; 
Lc_r = Lc_r/0.6*L;   % the interpolatin is based on Lc = 0.6; 
function z = func_H1(rsq,esq)

z = (2*esq + rsq)./(esq + rsq).^(3/2); 
z = z/(8*pi); 


function [angle] = cosine(l,m, c)



angle = (l'*c*m) / sqrt((l'*c*l)*(m'*c*m));



end
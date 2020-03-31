
function [ nor_par ] = norma_back( par, low, high )
%From parameters par in the range [low,upper], it gives the normalized
%parameters nor_par

if length(par)==1
    nor_par = (par-low)/(high-low);
else
    nor_par = (par-low)./(high-low);

end


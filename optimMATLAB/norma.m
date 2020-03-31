function [ par ] = norma( nor_par, low, high )
%From normalized parameter nor_par, it gives the parameter par in the range [low,upper]

if length(nor_par)==1
    par = nor_par*(high-low)+low;
else
    par = nor_par.*(high-low)+low;

end


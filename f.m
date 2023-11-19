function val = f(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%val = (1+5*pi*pi)*cos(pi*x).*cos(2*pi*y);
%val = cos(pi*x).*cos(2*pi*y) + 10*pi^2*cos(pi*x).*cos(2*pi*y) + pi^2*cos(2*pi*x).*sin(pi*x).*sin(4*pi*y) + 9*pi^2*cos(pi*x).*sin(2*pi*x).*cos(2*pi*y).*sin(2*pi*y);
%val = (1 + 5*pi^2) * sin(pi*x) .* sin(2*pi*y);
val = (1 + 5*pi^2) * cos(pi*x) .* cos(2*pi*y);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

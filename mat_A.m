function val = mat_A(x,y)%,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%val = [1 0;0 1];
%val = [1 0;0 2];
%val = [2+sin(2*pi*x) 0;0 4+sin(2*pi*x)];
val = [2+sin(2*pi*x) 0;0 4];
%val = sin(epsilon*2*pi*x)*sin(epsilon*2*pi*y) + 2;
%val = (2+sin(2*pi*x))*(4+sin(2*pi*y));
%val = [2+sin(2*pi*x/10) 0;0 2+sin(2*pi*x/10)];

%val = [2+sin(2*pi*x/100) 0;0 2+sin(2*pi*x/100)];
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

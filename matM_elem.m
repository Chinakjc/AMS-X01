function Mel = matM_elem(S1, S2, S3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matM_elem :
% calcul la matrices de masse elementaire en P1 lagrange
%
% SYNOPSIS Mel = matM_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Mel matrice de masse elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

%xs = [x1,x2,x3];
%ys = [y1,y2,y3];
% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

x_ij = zeros(3);
y_ij = zeros(3);

%for i=1:3
%    for j=1:3
%	    x_ij(i,j) = xs(i) - xs(j);
%		y_ij(i,j) = ys(i) - ys(j);


% calcul de la matrice de masse
% -----------------------------
%Mel = zeros(3,3);
%for i=1:3
%	for j=1:3
		% A COMPLETER
%		Mel(i,j) = 
%	end; % j
%end; % i

Mel = abs(D)/24.0*(eye(3)+ones(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre_per_cellule.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

%valeur de eta
eta = 2^(-100);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
%LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  tri = Numtri(l,:);
  S1=Coorneu(tri(1),:);
  S2=Coorneu(tri(2),:);
  S3=Coorneu(tri(3),:);
  % calcul des matrices elementaires du triangle l 
  
   %Kel=matKvar_elem(S1, S2, S3);
   Kel=matK_elem(S1, S2, S3);
           
   Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemblage de la matrice globale et du second membre
  % A COMPLETER
   for i=1:3
      I = tri(i);
      for j=1:3
          J = tri(j);
          MM(I,J) += Mel(i,j);
          KK(I,J) += Kel(i,j);
      end
  end      

end % for l

% Calcul du second membre L
%(x1,y1),...,(xn,yn)] est la liste des cornneu, alors XX = [x1,...,xn] et YY = [y1,...,yn]
XX = Coorneu(:,1);
YY = Coorneu(:,2);

%e1 = grad XX et e2 = grad YY, LLX = -K*XX et LLY = -K*YY
LLX = -KK*XX;
LLY = -KK*YY;

% Projection sur l espace V_p
% ———————————————————
% matrice de projection 
PP = sparse(Nbpt,Nbpt);
Nbpt_interieur = 0;
for indice = 1:Nbpt
    if Refneu(indice) == 0
        PP(indice,indice) = 1;
        Nbpt_interieur += 1;
    end
end %for indice
Nbpt_bord = Nbpt - Nbpt_interieur;
PP(1,1) = 1;
PP(1,2) = 1;
PP(1,3) = 1;
PP(1,4) = 1;
Nbpt_bord_1 = 0;
for indice = 5:Nbpt_bord
    if Coorneu(indice,2) == 0
       Nbpt_bord_1 += 1;
    else 
        break;
    end %if
end   %for indice

Nbpt_bord_2 = Nbpt_bord/2 - Nbpt_bord_1 - 2;
for indice = 5: 4 + Nbpt_bord_1 
    PP(indice,indice) = 1;
    PP(5 + Nbpt_bord_1 - indice + 4,indice + Nbpt_bord_1 + Nbpt_bord_2) = 1;

end %for indice

for indice = 5+Nbpt_bord_1:4 + Nbpt_bord_1+Nbpt_bord_2 
    PP(indice,indice) = 1;
    PP(5 + Nbpt_bord_1 + Nbpt_bord_2 + Nbpt_bord_1 - indice + 4,indice + Nbpt_bord_1 + Nbpt_bord_2) = 1;
end


AA = eta*MM+KK;
AAp = PP*AA*PP';
LLXp = PP*LLX;
LLYp = PP*LLY;

% inversion
% ----------
WXp = AAp\LLXp;
WYp = AAp\LLYp;

% Expression de la solution dans toute la base
% ———————
W1 = PP'*WXp;
W2 = PP'*WYp;


% visualisation
% -------------
affiche(W1, Numtri, Coorneu, sprintf('Cellule i=1 - %s', nom_maillage));
affiche(W2, Numtri, Coorneu, sprintf('Cellule i=2 - %s', nom_maillage));

%tenseur homogeneise
Aeff = zeros(2,2);

Aeff(1,1) = (XX + W1)' * (KK * (XX + W1));
Aeff(1,2) = (XX + W1)' * (KK * (YY + W2));
Aeff(2,1) = (YY + W2)' * (KK * (XX + W1));
Aeff(2,2) = (YY + W2)' * (KK * (YY + W2));

Aeff 

validation = 'non';
% validation
% ----------
if strcmp(validation,'oui')
W1_exact = zeros(size(Coorneu));%cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
W2_exact = zeros(size(Coorneu));
%affiche(UU_exact, Numtri, Coorneu, sprintf('Periodique -exact %s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
errX = W1 - W1_exact;
errY = W2 - W2_exact;
errL2 = sqrt(errX' * (MM * errX) + errY' * (MM * errY));
% Calcul de l erreur H1
errH1 = sqrt(errL2^2 + errX' * (KK*errX) + errY' * (KK*errY));
% A COMPLETER
printf("erreur L2 = %f et erreur H1 = %f", errL2,errH1);
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
nom_maillage = 'geomCarre_per.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

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
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

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


AA = MM+KK;
AAp = PP*AA*PP';
LLp = PP*LL;

% inversion
% ----------
UUp = AAp\LLp;

% Expression de la solution dans toute la base
% ———————
UU = PP'*UUp;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = cos(pi*Coorneu(:,1)).*cos(2*pi*Coorneu(:,2));
affiche(UU_exact, Numtri, Coorneu, sprintf('Periodique -exact %s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
errU = UU_exact - UU;
errL2 = sqrt(errU' * (MM * errU));
% Calcul de l erreur H1
errH1 = sqrt(errL2^2 + errU' * (KK*errU));
% A COMPLETER
printf("erreur L2 = %f et erreur H1 = %f", errL2,errH1);
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


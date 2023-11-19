% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'geomCarre.msh';
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
  
   Kel=matK_elem(S1, S2, S3);
           
   Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
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

% Projection sur l espace V_0
% ———————————————————
% matrice de projection 
Point_interieur = find(Refneu == 0);
Nbpt_interieur = size(Point_interieur,1);
PP = sparse(Nbpt_interieur,Nbpt);
for indice = 1:Nbpt_interieur
    PP(indice,Point_interieur(indice)) = 1;
end
AA = MM+KK;
AA0 = PP*AA*PP';
LL0 = PP*LL;

% inversion
% ----------
UU0 = AA0\LL0;

% Expression de la solution dans toute la base
% ———————
UU = PP'*UU0;

% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
%UU_exact = sin(pi*Coorneu(:,1)).*sin(2*pi*Coorneu(:,2));
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
affiche(UU_exact, Numtri, Coorneu, sprintf('Dirichlet -exact %s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
errU = UU_exact - UU;
errU0 = PP*errU;
errL2 = sqrt(errU' * (MM * errU));
%errL2 = sqrt(errU0' * (PP*MM*PP' * errU0));
% Calcul de l erreur H1
errH1 = sqrt(errL2^2 + errU' * (KK*errU));
%errH1 = sqrt(errU0' * AA0 * errU0);
% A COMPLETER
printf("erreur L2 = %f et erreur H1 = %f", errL2,errH1);
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


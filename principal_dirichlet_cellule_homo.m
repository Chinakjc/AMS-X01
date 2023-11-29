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

func_u_exact = @(x,y) sin(2*pi*x).*sin(2*pi*y);

%Exo2 Question11 
%Exo 3 Question 1,2,3,4
func_A_list = cell(1, 5); 

func_A_list{1} = @(x,y) 1;
func_A_list{2} = @(x,y) [1, 0; 0, 2];
func_A_list{3} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4];
func_A_list{4} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4 + sin(2*pi*x)];
func_A_list{5} = @(x,y) (2 + sin(2*pi*x)).*(4 + sin(2*pi*y));

func_f_list = cell(1, 5);

func_f_list{1} = @(x,y) 8*pi^2*sin(2*pi*x).*sin(2*pi*y);
func_f_list{2} = @(x,y) 12*pi^2*sin(2*pi*x).*sin(2*pi*y);
func_f_list{3} = @(x,y) -4*pi^2*cos(2*pi*x).*cos(2*pi*x).*sin(2*pi*y) + 4*pi^2*(sin(2*pi*x)+ 2).*sin(2*pi*x).*sin(2*pi*y) + 16*pi^2*sin(2*pi*x).*sin(2*pi*y);
func_f_list{4} = @(x,y) -4*pi^2*cos(2*pi*x).*cos(2*pi*x).*sin(2*pi*y) + 4*pi^2*(sin(2*pi*x)+ 4).*sin(2*pi*x).*sin(2*pi*y) + 4*pi^2*(sin(2*pi*x)+ 2).*sin(2*pi*x).*sin(2*pi*y);
func_f_list{5} = @(x,y) -4*pi^2*(sin(2*pi*x)+ 2).*cos(2*pi*y).*cos(2*pi*y).*sin(2*pi*x) - 4*pi^2*(sin(2*pi*y)+ 4).*cos(2*pi*x).*cos(2*pi*x).*sin(2*pi*y) + 8*pi^2*(sin(2*pi*x)+ 2).*(sin(2*pi*y)+ 4).*sin(2*pi*x).*sin(2*pi*y);

index_problem = 5; % on choisit le probleme a resoudre parmi les 5

validation = 'oui';



func_A = func_A_list{index_problem};
func_f = func_f_list{index_problem};

% calcul du du tenseur homogénéisé
% --------------------------------

%valeur de eta
eta = 2^(-8);

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'validation/geomCarre_per_cellule.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite

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
   Kel=matK_elem_dynamique(func_A,S1, S2, S3);
           
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
% fin du calcul du tenseur homogénéisé
% ------------

% calcul de la solution homogénéisée
% ----------------------------------

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'validation/geomCarre4.msh';
fprintf('%s\n', nom_maillage);
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
func_Aeff = @(x,y) Aeff;

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
  
   Kel=matK_elem_dynamique(func_A,S1, S2, S3);
           
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
FF = func_f(Coorneu(:,1),Coorneu(:,2));
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
AA = KK;
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
affiche(UU, Numtri, Coorneu, sprintf('Dirichlet omogénéisé- %s', nom_maillage));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
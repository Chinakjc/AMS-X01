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

%Exo 2 Question 7,8
func_w1_exact = @(x,y) zeros(size(x));
func_w2_exact = @(x,y) zeros(size(x));

validation = 'non';%'oui';

%valeur de eta
eta = 2^(-8);

%Exo 2 Question 7,8,9,10
if strcmp(validation,'oui')
%Exo 2 Question 7,8
func_A_list = cell(1, 2); 

func_A_list{1} = @(x,y) 1;
func_A_list{2} = @(x,y) [1, 0; 0, 2];
else
%Exo 2 Question 9,10
func_A_list = cell(1, 5); 

func_A_list{1} = @(x,y) 1;
func_A_list{2} = @(x,y) [1, 0; 0, 2];
func_A_list{3} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4];
func_A_list{4} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4 + sin(2*pi*x)];
func_A_list{5} = @(x,y) (2 + sin(2*pi*x)).*(4 + sin(2*pi*y));

end


% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'validation/geomCarre_per_cellule.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

for index_func = 1:length(func_A_list)

func_A = func_A_list{index_func};

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

fprintf('cas %d\n', index_func - 1);
Aeff 


% validation
% ----------
if strcmp(validation,'oui')

W1_exact = func_w1_exact(Coorneu(:,1), Coorneu(:,2));
W2_exact = func_w2_exact(Coorneu(:,1), Coorneu(:,2));
%affiche(UU_exact, Numtri, Coorneu, sprintf('Periodique -exact %s', nom_maillage));
% Calcul de l erreur L2
% A COMPLETER
errX = W1 - W1_exact;
errY = W2 - W2_exact;
errL2 = sqrt(errX' * (MM * errX) + errY' * (MM * errY));
% Calcul de l erreur H1
errH1 = sqrt(errL2^2 + errX' * (KK*errX) + errY' * (KK*errY));
% A COMPLETER
fprintf("erreur L2 = %f et erreur H1 = %f\n", errL2,errH1);
% attention de bien changer le terme source (dans FF)
end

end %for index_func

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

func_u_exact = @(x,y) sin(pi*x).*sin(pi*y);

%Exo 2 Question 2,3,4
func_A_list = cell(1, 5); 

func_A_list{1} = @(x,y) 1;
func_A_list{2} = @(x,y) [1, 0; 0, 2];
func_A_list{3} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4];
func_A_list{4} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4 + sin(2*pi*x)];
func_A_list{5} = @(x,y) (2 + sin(2*pi*x)).*(4 + sin(2*pi*y));

func_f_list = cell(1, 5);

func_f_list{1} = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
func_f_list{2} = @(x,y) 3*pi^2*sin(pi*x).*sin(pi*y);
func_f_list{3} = @(x,y) -2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y) + pi^2*(sin(2*pi*x)+ 2).*sin(pi*x).*sin(pi*y) + 4*pi^2*sin(pi*x).*sin(pi*y);
func_f_list{4} = @(x,y) -2*pi^2*cos(2*pi*x).*cos(pi*x).*sin(pi*y) + pi^2*(sin(2*pi*x)+ 4).*sin(pi*x).*sin(pi*y) + pi^2*(sin(2*pi*x)+ 2).*sin(pi*x).*sin(pi*y);
func_f_list{5} = @(x,y) -2*pi^2*(sin(2*pi*x)+ 2).*cos(2*pi*y).*cos(pi*y).*sin(pi*x) - 2*pi^2*(sin(2*pi*y)+ 4).*cos(2*pi*x).*cos(pi*x).*sin(pi*y) + 2*pi^2*(sin(2*pi*x)+ 2).*(sin(2*pi*y)+ 4).*sin(pi*x).*sin(pi*y);


validation = 'oui';

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'validation/geomCarre5.msh';
fprintf('%s\n', nom_maillage);
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------



UU_exact = [];
if strcmp(validation,'oui')
UU_exact = func_u_exact(Coorneu(:,1),Coorneu(:,2));
affiche(UU_exact, Numtri, Coorneu, sprintf('Dirichlet -exact %s', nom_maillage));
end


for index_func = 1:length(func_A_list)

func_A = func_A_list{index_func};
func_f = func_f_list{index_func};

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


% validation
% ----------
if strcmp(validation,'oui')

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
fprintf("erreur L2 = %f et erreur H1 = %f\n", errL2,errH1);
% attention de bien changer le terme source (dans FF)



end

end % for index_func


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
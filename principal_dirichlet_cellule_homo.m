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
func_uxx_exact = @(x,y) -4*pi^2*sin(2*pi*x).*sin(2*pi*y);
func_uyy_exact = @(x,y) -4*pi^2*sin(2*pi*x).*sin(2*pi*y);
func_uxy_exact = @(x,y) 4*pi^2*cos(2*pi*x).*cos(2*pi*y);

%Exo2 Question11 
%Exo 3 Question 1,2,3,4
func_A_list = cell(1, 5); 

func_A_list{1} = @(x,y) 1;
func_A_list{2} = @(x,y) [1, 0; 0, 2];
func_A_list{3} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4];
func_A_list{4} = @(x,y) [2 + sin(2*pi*x), 0; 0, 4 + sin(2*pi*x)];
func_A_list{5} = @(x,y) (2 + sin(2*pi*x)).*(4 + sin(2*pi*y));

%func_f_list = cell(1, 5);

%func_f_list{1} = @(x,y) 8*pi^2*sin(2*pi*x).*sin(2*pi*y);
%func_f_list{2} = @(x,y) 12*pi^2*sin(2*pi*x).*sin(2*pi*y);
%func_f_list{3} = @(x,y) -4*pi^2*cos(2*pi*x).*cos(2*pi*x).*sin(2*pi*y) + 4*pi^2*(sin(2*pi*x)+ 2).*sin(2*pi*x).*sin(2*pi*y) + 16*pi^2*sin(2*pi*x).*sin(2*pi*y);
%func_f_list{4} = @(x,y) -4*pi^2*cos(2*pi*x).*cos(2*pi*x).*sin(2*pi*y) + 4*pi^2*(sin(2*pi*x)+ 4).*sin(2*pi*x).*sin(2*pi*y) + 4*pi^2*(sin(2*pi*x)+ 2).*sin(2*pi*x).*sin(2*pi*y);
%func_f_list{5} = @(x,y) -4*pi^2*(sin(2*pi*x)+ 2).*cos(2*pi*y).*cos(2*pi*y).*sin(2*pi*x) - 4*pi^2*(sin(2*pi*y)+ 4).*cos(2*pi*x).*cos(2*pi*x).*sin(2*pi*y) + 8*pi^2*(sin(2*pi*x)+ 2).*(sin(2*pi*y)+ 4).*sin(2*pi*x).*sin(2*pi*y);

index_problem = 5; % on choisit le probleme a resoudre parmi les 5

validation = 'oui';
validation_osc = 'oui'; % Exo 3 Question 1,2,3,4


func_A = func_A_list{index_problem};
%func_f = func_f_list{index_problem};

% calcul du du tenseur homogénéisé
% --------------------------------
fprintf('Calcul du tenseur homogénéisé\n');
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
fprintf('fin du calcul du tenseur homogénéisé\n');
% fin du calcul du tenseur homogénéisé
% ------------

% calcul de la solution homogénéisée
% ----------------------------------
fprintf('Calcul de la solution homogénéisée\n');

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'validation/geomCarre5.msh';
fprintf('%s\n', nom_maillage);
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
func_Aeff = @(x,y) Aeff;
func_f = @(x,y) -(Aeff(1,1)*func_uxx_exact(x,y) + Aeff(1,2)*func_uxy_exact(x,y) + Aeff(2,1)*func_uxy_exact(x,y) + Aeff(2,2)*func_uyy_exact(x,y));

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
  
   Kel=matK_elem_dynamique(func_Aeff,S1, S2, S3);
           
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

fprintf('fin du calcul de la solution homogénéisée\n');
% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Dirichlet homogénéisé- %s', nom_maillage));

% validation
% ----------

if strcmp(validation,'oui')
UU_exact = func_u_exact(Coorneu(:,1),Coorneu(:,2));
affiche(UU_exact, Numtri, Coorneu, sprintf('Dirichlet homogénéisé -exact %s', nom_maillage));
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
end

if strcmp(validation_osc,'oui')

% calcul de la solution ocsillante
% ---------------------------------

%Exo 3 Question 1,2,3,4
epsilons = 8.^[0:1:10];

errL2_list = [];
errH1_list = [];

func_A_eps_list = cell(1, size(epsilons,2));
for index_eps = 1:size(epsilons,2)
epsilon = epsilons(index_eps);
func_A_eps_list{index_eps} = @(x,y) func_A(x/epsilon,y/epsilon);
end


fprintf('Calcul de la solution ocsillante\n');

% lecture du maillage et affichage
% ---------------------------------
nom_maillage = 'validation/geomCarre5.msh';
fprintf('%s\n', nom_maillage);
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);


for index_eps = 1:size(epsilons,2)

epsilon = epsilons(index_eps);
fprintf('epsilon = %d\n', epsilon);

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre
func_A_osc = func_A_eps_list{index_eps};

KK_std = sparse(Nbpt,Nbpt); % matrice de rigidite standard

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
  
   Kel=matK_elem_dynamique(func_A_osc,S1, S2, S3);
   Kel_std = matK_elem_old(S1, S2, S3);
           
   Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  % A COMPLETER
   for i=1:3
      I = tri(i);
      for j=1:3
          J = tri(j);
          MM(I,J) += Mel(i,j);
          KK(I,J) += Kel(i,j);
          KK_std(I,J) += Kel_std(i,j);
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
UU_osc = PP'*UU0;

fprintf('fin du calcul de la solution ocsillante\n');
% visualisation
% -------------
affiche(UU_osc, Numtri, Coorneu, sprintf('Dirichlet ocsillant- %s avec epsilon = %d', nom_maillage,epsilon));

% validation
errU_osc = UU - UU_osc;
errL2_osc = sqrt(errU_osc' * (MM * errU_osc));
errH1_osc = sqrt(errL2_osc^2 + errU_osc' * (KK_std*errU_osc));
fprintf("erreur L2 = %f et erreur H1 = %f\n", errL2_osc,errH1_osc);
errL2_list = [errL2_list, errL2_osc];
errH1_list = [errH1_list, errH1_osc];

end %for index_eps

% fin du calcul de la solution ocsillante
% ---------------------------------------


% affichage de l'erreur en fonction de epsilon
figure();

% affichage de l'erreur en norme L2
scatter(-log2(epsilons), log2(errL2_list), 'r', 'DisplayName', 'erreur L2');
hold on;
% affichage de l'erreur en norme H1
scatter(-log2(epsilons), log2(errH1_list), 'b', 'DisplayName', 'erreur H1');

% asymptote L2
x_last_two = log2(epsilons(end-1:end));
y_last_two_L2 = log2(errL2_list(end-1:end));

slope_L2 = (y_last_two_L2(2) - y_last_two_L2(1)) / (x_last_two(2) - x_last_two(1));
intercept_L2 = y_last_two_L2(2) - slope_L2 * x_last_two(2);

x_asymptote_L2 = linspace(min(-log2(epsilons)), max(-log2(epsilons)), 128);
y_asymptote_L2 = -slope_L2 * x_asymptote_L2 + intercept_L2;

plot(x_asymptote_L2, y_asymptote_L2, '--r', 'DisplayName', 'Asymptote L2');

% asymptote H1
y_last_two_H1 = log2(errH1_list(end-1:end));

slope_H1 = (y_last_two_H1(2) - y_last_two_H1(1)) / (x_last_two(2) - x_last_two(1));
intercept_H1 = y_last_two_H1(2) - slope_H1 * x_last_two(2);

y_asymptote_H1 = -slope_H1 * x_asymptote_L2 + intercept_H1;

plot(x_asymptote_L2, y_asymptote_H1, '--b', 'DisplayName', 'Asymptote H1');

hold off;

legend();
title('Erreur en fonction de 1/epsilon en échelle log2');

% changement des echelles des axes
set(gca, 'xticklabel', arrayfun(@(x) ['2^{' num2str(x) '}'], get(gca, 'xtick'), 'UniformOutput', false));
set(gca, 'yticklabel', arrayfun(@(y) ['2^{' num2str(y) '}'], get(gca, 'ytick'), 'UniformOutput', false));

% changement de la taille de la figure
set(gcf, 'Position', [100, 100, 1280, 1024]);

% ordre de convergence
fprintf('Ordre de convergence pour L2: %f\n', -slope_L2);
fprintf('Ordre de convergence pour H1: %f\n', -slope_H1);


end %if strcmp(validation_osc,'oui')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
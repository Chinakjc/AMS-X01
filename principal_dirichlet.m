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

func_u_exact = @(x,y) sin(2*pi*x).*sin(2*pi*y);%Question 16
func_f = @(x,y) (1+8*pi*pi)*sin(2*pi*x).*sin(2*pi*y); %Question 16
func_A = @(x,y) 1; %Question 16


validation = 'oui'; %Question 16


errsL2 = [];
errsH1 = [];


if strcmp(validation,'oui')
maillages = {"validation/geomCarre2.msh", "validation/geomCarre3.msh", "validation/geomCarre4.msh", "validation/geomCarre5.msh"};
hs = [1/4, 1/8, 1/16, 1/32];
else
maillages = {'validation/geomCarre5.msh'};
hs = [1/32];
end


for nom_maillage = maillages
% lecture du maillage et affichage
% ---------------------------------
nom_maillage = char(nom_maillage);
fprintf('%s\n', nom_maillage);
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
UU_exact = func_u_exact(Coorneu(:,1),Coorneu(:,2));
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
fprintf("erreur L2 = %f et erreur H1 = %f\n", errL2,errH1);
% attention de bien changer le terme source (dans FF)

errsL2 = [errsL2, errL2];
errsH1 = [errsH1, errH1];

end

end % for nom_maillage


if strcmp(validation,'oui')

figure();

% affichage de l'erreur en norme L2
scatter(-log2(hs), log2(errsL2), 'r', 'DisplayName', 'erreur L2');
hold on;
% affichage de l'erreur en norme H1
scatter(-log2(hs), log2(errsH1), 'b', 'DisplayName', 'erreur H1');

% asymptote L2
x_last_two = log2(hs(end-1:end));
y_last_two_L2 = log2(errsL2(end-1:end));

slope_L2 = (y_last_two_L2(2) - y_last_two_L2(1)) / (x_last_two(2) - x_last_two(1));
intercept_L2 = y_last_two_L2(2) - slope_L2 * x_last_two(2);

x_asymptote_L2 = linspace(min(-log2(hs)), max(-log2(hs)), 128);
y_asymptote_L2 = -slope_L2 * x_asymptote_L2 + intercept_L2;

plot(x_asymptote_L2, y_asymptote_L2, '--r', 'DisplayName', 'Asymptote L2');

% asymptote H1
y_last_two_H1 = log2(errsH1(end-1:end));

slope_H1 = (y_last_two_H1(2) - y_last_two_H1(1)) / (x_last_two(2) - x_last_two(1));
intercept_H1 = y_last_two_H1(2) - slope_H1 * x_last_two(2);

y_asymptote_H1 = -slope_H1 * x_asymptote_L2 + intercept_H1;

plot(x_asymptote_L2, y_asymptote_H1, '--b', 'DisplayName', 'Asymptote H1');

hold off;

legend();
title('Erreur en fonction de 1/h en échelle log2');

% changement des echelles des axes
set(gca, 'xticklabel', arrayfun(@(x) ['2^{' num2str(x) '}'], get(gca, 'xtick'), 'UniformOutput', false));
set(gca, 'yticklabel', arrayfun(@(y) ['2^{' num2str(y) '}'], get(gca, 'ytick'), 'UniformOutput', false));

% changement de la taille de la figure
set(gcf, 'Position', [100, 100, 1280, 1024]);

% ordre de convergence
fprintf('Ordre de convergence pour L2: %f\n', slope_L2);
fprintf('Ordre de convergence pour H1: %f\n', slope_H1);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


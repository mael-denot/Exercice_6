% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%
clear; clc; close all;
%% Parametres %%
%%%%%%%%%%%%%%%%

%Chemin d'acces au code compile
repertoire = './'; % './' on Linux, '' on Windows
executable = 'Exercice6_2023_student.exe'; % Nom de l'executable


input = 'configuration.in.example';
N1       = 10:2:20;
N2       = N1;
nsimul   = numel(N1);
paramstr = 'N'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param    = [N1;N2]; % Valeurs du parametre a scanner


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul);

for ii = 1:nsimul
    % Variant to scan N1 and N2 together:
    filename  = [paramstr, '=', num2str(param(1,ii))];
    output{ii} = [filename, '.out'];
    eval(sprintf('!%s%s %s %s=%.15g %s=%.15g output=%s', repertoire, executable, input, [paramstr,'1'], param(1,ii), [paramstr,'2'], param(2,ii), output{ii}));
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

for ii = 1:nsimul
    file_phi   = [output{ii},'_phi.out'];
    data_phi   = load(file_phi);
end


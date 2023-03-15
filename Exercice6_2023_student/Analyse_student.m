%% Chargement des résultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier_phi = 'output_phi.out';
fichier_E   = 'output_E.out';
fichier_D   = 'output_D.out';
data = load(fichier_phi);
r=data(:,1);
phi=data(:,2);

data = load(fichier_E);
rmid=data(:,1);
E=data(:,2);

data = load(fichier_D);
D=data(:,2);

%% Figures %%
%%%%%%%%%%%%%
figure
h = plotyy(r,phi,rmid,E);
xlabel('r [m]')
ylabel(h(1),'\phi [V]')
ylabel(h(2),'E_{r} [V/m]')
title('Potentiel electrique \phi=\phi(r) et champ electrique E_{r}=E_{r}(r)')

grid

figure
h = plotyy(rmid,D,rmid,E);
xlabel('r [m]')
ylabel(h(1),'D_{r}/ \epsilon_{0} [V/m]')
ylabel(h(2),'E_{r} [V/m]')
title('champ de déplacement D_{r}=D_{r}(r) et champ electrique E_{r}=E_{r}(r)')

grid


%% Verification of continuity equation
%%%%%%%%%%%%%%%

%% TO DO: Compute div_D with forward finite differences
ii  = 1:(length(rmid)-1);
div_D   = 0*ii;

%% TO D0: Compute the free charge profile /eps_0 with an handle func
%%        (See MATLAB documentation for handle func)         
rho = @(x) 0*x;

figure('Name','Equation de continuité')
plot(rmid(ii),div_D,'b.')
xlabel('r [m]')
ylabel('div(D)[C/m^{3}]')
grid
hold on
plot(r,rho(r),'r')
legend('div(D)','\rho(r)')
hold off


%% Code to make figures 7a,b,c,d


clear all
Tol = 0.000000000001;

% Possible Shears.
ROneSize = 1.1;
Phi1 = [1,1,1,1,1];

BAll = 0.001:0.0001:0.2;
for i=1:length(BAll)
    B=BAll(i);
Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
R = [ROneSize, Rtmp];
B1All = linspace(1.0001,1.5001,3000);
[Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All);

[MaxCp,MaxCPInd] = max(Cp);

KAll(i,:)=kk;

CtAll(i,:) = kk.*A2.^2;
CpAll(i,:) = Cp;
end

X = repmat(BAll',1,length(KAll(1,:)));
Y = KAll;
Z = CpAll;
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
[C,h] = contourf(X,Y,Z,10);
h.LevelList=round(h.LevelList,2);  %rounds levels to 3rd decimal places
%clabel(C,h,'FontSize',15)
ylim([0 16])
xlim([0.01 0.2])
c = colorbar;
ylabel(c,'$C_{P}^{*}$','Interpreter','latex');
caxis([0 0.85])
xlabel('$B$','Interpreter','latex')
ylabel('$k$','Interpreter','latex')
[~,tmp] = max(Z,[],2);
hold on
for i=1:length(BAll)
   Kmax(i) = Y(i,tmp(i)); 
end
plot(BAll,Kmax,'r','LineWidth',2)
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig7a.eps','Resolution',1200)
hold off

gamma = 2 ;
lamb = BAll/4;
cf0 = 0.002;
lmbda_cf0 = lamb./cf0;



zeta = 5;
a=(CtAll.*lmbda_cf0' +1);
b=zeta;
c= -1-zeta;
pos = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
CppAll=(pos.^3).*CpAll;
Z = CppAll;
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
[C,h] = contourf(X,Y,Z,10);
h.LevelList=round(h.LevelList,2);  %rounds levels to 3rd decimal places

ylim([0 16])
xlim([0.01 0.2])
c = colorbar;
ylabel(c,'$C_{P}^{G}$','Interpreter','latex');
caxis([0 0.6])
xlabel('$B$','Interpreter','latex')
ylabel('$k$','Interpreter','latex')
[~,tmp] = max(Z,[],2);
hold on
for i=1:length(BAll)
   Kmax(i) = Y(i,tmp(i)); 
end
plot(BAll,Kmax,'r','LineWidth',2)
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig7b.eps','Resolution',1200)
hold off

zeta = 15;
a=(CtAll.*lmbda_cf0' +1);
b=zeta;
c= -1-zeta;
pos = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
CppAll=(pos.^3).*CpAll;
Z = CppAll;
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
[C,h] = contourf(X,Y,Z,10);
h.LevelList=round(h.LevelList,2);  %rounds levels to 3rd decimal places

ylim([0 16])
xlim([0.01 0.2])
c = colorbar;
ylabel(c,'$C_{P}^{G}$','Interpreter','latex');
caxis([0 0.6])
xlabel('$B$','Interpreter','latex')
ylabel('$k$','Interpreter','latex')
[~,tmp] = max(Z,[],2);
hold on
for i=1:length(BAll)
   Kmax(i) = Y(i,tmp(i)); 
end
plot(BAll,Kmax,'r','LineWidth',2)
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig7c.eps','Resolution',1200)
hold off


zeta = 25;
a=(CtAll.*lmbda_cf0' +1);
b=zeta;
c= -1-zeta;
pos = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
CppAll=(pos.^3).*CpAll;

Z = CppAll;

figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
[C,h] = contourf(X,Y,Z,10);
h.LevelList=round(h.LevelList,2);  %rounds levels to 3rd decimal places
ylim([0 16])
xlim([0.01 0.2])
c = colorbar;
ylabel(c,'$C_{P}^{G}$','Interpreter','latex');
caxis([0 0.6])

xlabel('$B$','Interpreter','latex')
ylabel('$k$','Interpreter','latex')
[~,tmp] = max(Z,[],2);
hold on
for i=1:length(BAll)
   Kmax(i) = Y(i,tmp(i)); 
end
plot(BAll,Kmax,'r','LineWidth',2)
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig7d.eps','Resolution',1200)
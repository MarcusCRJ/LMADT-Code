%% Code for figure 8

clear all


Tol = 0.000000000001;

%figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])

% Possible Shears.
ROneSize = 1.1;
VelAv = 1;


V = VelAv;
Phi1 = [1,V,V,V,V];

BAll = 0.001:0.001:0.2;
CpAllMax = zeros(1,length(BAll));
CtAllMax = CpAllMax;
for i=1:length(BAll)
    B=BAll(i);
Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
R = [ROneSize, Rtmp];
B1All = 1.001:0.001:10;
[Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All);

[MaxCp,MaxCPInd] = max(Cp);

CtAllMax(i) = kk(MaxCPInd)*A2(MaxCPInd)^2;
CpAllMax(i) = MaxCp;
end
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
hold on
plot(BAll,CpAllMax,'k')
gamma = 2 ;
%#lamb = 50/(75*20*100) # D / 20D * (Width/50)
%width = df.iloc[j,2]
%lamb = 100./(20*100.*(1./BAll));
%lamb = BAll/20;
lamb = BAll/4;
cf0 = 0.002;
lmbda_cf0 = lamb./cf0;

zeta = 5;
a=(CtAllMax.*lmbda_cf0 +1);
b=zeta;
c= -1-zeta;
pos = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
CppAll=(pos.^3).*CpAllMax;
plot(BAll,CppAll,'k-.')


zeta = 15;
a=(CtAllMax.*lmbda_cf0 +1);
b=zeta;
c= -1-zeta;
pos = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
CppAll=(pos.^3).*CpAllMax;
plot(BAll,CppAll,'k-*')


zeta = 25;
a=(CtAllMax.*lmbda_cf0 +1);
b=zeta;
c= -1-zeta;
pos = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
CppAll=(pos.^3).*CpAllMax;
plot(BAll,CppAll,'k--')

xlim([0 0.2])

legend('$C_{P\:max}^{*}$','$C_{P}^{G}$ at $k=k_{opt}^{*}$, $\zeta = 5$','$C_{P}^{G}$ at $k=k_{opt}^{*}$, $\zeta = 15$','$C_{P}^{G}$ at $k=k_{opt}^{*}$, $\zeta = 25$','Uniform Flow','Interpreter','latex')
xlabel('$B$','Interpreter','latex')
ylabel('Power Coefficient','Interpreter','latex')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig8.eps','Resolution',1200)

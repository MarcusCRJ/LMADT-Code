%% Code for figure 5
clear all
tic
Tol = 0.000000000001;
B=0.2;

RAlls = 1;
KAll = 0.1:0.1:16;
VelAv = [0.501,1,1.5,2,5,10];
ResAll = linspace(1,400,400);
CpAll = zeros(length(ResAll),length(VelAv));


for ii=1:length(ResAll)
    for k=1:length(VelAv)
        V = VelAv(k);
        Phi1 = [1,linspace(1,2*V-1,ResAll(ii))];
        CpTmp = zeros(length(KAll),1);
        ROneSize = RAlls;
        Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
        R = [ROneSize, Rtmp];
        B1All = 1.05:0.05:10;
        [Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All);
        C = max(Cp(:));
        CpAll(ii,k) = C;
    end
end
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
plot(ResAll,CpAll)
legend('$\phi=0.5$','$\phi=1$','$\phi=1.5$','$\phi=2$','$\phi=5$','$\phi=10$','Interpreter','latex')
toc
xlabel('$N$','Interpreter','latex')
ylabel('$C_{P\:\textit{max}}$','Interpreter','latex')
ylim([0 6])
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig5.eps','Resolution',1200)
%% Make figure 4 plots

clear all
tic
Tol = 0.000000000001;
B=0.2;
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])

% Possible Shears.
RAlls = linspace(1,4.5,81);
KAll = 0.1:0.1:16;
VelAv = linspace(0.51,10,151);

%VelAv = linspace(0.1,1300,81);

EffBlock1 = zeros(length(RAlls),length(VelAv));
EffBlock2 = zeros(length(RAlls),length(VelAv));
CpAll = zeros(length(KAll),length(RAlls),length(VelAv));

for k=1:length(VelAv)
    C = zeros(length(RAlls),1);
    V = VelAv(k);
    Phi1 = [1,V,V,V,V];
    CpTmp = zeros(length(KAll),1);
    for j=1:length(RAlls)
        ROneSize = RAlls(j);
        Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
        R = [ROneSize, Rtmp];
        B1All = 1.01:0.01:10;
        [Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All,0);
        C(j) = max(Cp(:));
    end
    EffBlock1(:,k) = 1 - sqrt(16./(27*C));
end
toc

levels = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
levelsCols = ["r";"#A2142F";"#D95319";"b";"#0072BD";"c";"g";"#77AC30";"k"];
for jj=1:length(levels)
    hold on
    [C,h] = contour(VelAv,RAlls,EffBlock1,[levels(jj) levels(jj)]);
    
    h.LineColor = levelsCols(jj);
    clabel(C,h,'FontSize',15)
    hold off
end
xlabel('$\phi$','Interpreter','latex')
ylabel('$r$','Interpreter','latex')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Figure4a.eps','Resolution',1200)
hold off

tic
for k=1:length(VelAv)
    C = zeros(length(RAlls),1);
    V = VelAv(k);
    Phi1 = [1,linspace(1,2*V-1,200)];
    
    CpTmp = zeros(length(KAll),1);
    for j=1:length(RAlls)
        ROneSize = RAlls(j);
        Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
        R = [ROneSize, Rtmp];
        B1All = 1.05:0.05:10;
        [Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All);
        C(j) = max(Cp(:));
    end
    EffBlock2(:,k) = 1 - sqrt(16./(27*C));
end
toc
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])

for jj=1:length(levels)
    hold on
    [C,h] = contour(VelAv,RAlls,EffBlock2,[levels(jj) levels(jj)]);
    
    h.LineColor = levelsCols(jj);
    clabel(C,h,'FontSize',15)
    hold off
end
%clabel(C,h,'FontSize',15)
xlabel('$\phi$','Interpreter','latex')
ylabel('$r$','Interpreter','latex')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Figure4b.eps','Resolution',1200)


B = 0.1;

hold off
tic
for k=1:length(VelAv)
    C = zeros(length(RAlls),1);
    V = VelAv(k);
    Phi1 = [1,linspace(1,2*V-1,200)];
    
    CpTmp = zeros(length(KAll),1);
    for j=1:length(RAlls)
        ROneSize = RAlls(j);
        Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
        R = [ROneSize, Rtmp];
        B1All = 1.05:0.05:10;
        [Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All);
        C(j) = max(Cp(:));
    end
    EffBlock2(:,k) = 1 - sqrt(16./(27*C));
end
toc
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])

for jj=1:length(levels)
    hold on
    [C,h] = contour(VelAv,RAlls,EffBlock2,[levels(jj) levels(jj)]);
    
    h.LineColor = levelsCols(jj);
    clabel(C,h,'FontSize',15)
    hold off
end
xlabel('$\phi$','Interpreter','latex')
ylabel('$r$','Interpreter','latex')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Figure4c.eps','Resolution',1200)

hold off

for k=1:length(VelAv)
    C = zeros(length(RAlls),1);
    V = VelAv(k);
    Phi1 = [1,V,V,V,V];
    CpTmp = zeros(length(KAll),1);
    for j=1:length(RAlls)
        ROneSize = RAlls(j);
        Rtmp = (1/B - ROneSize)*ones(1,length(Phi1)-1)/(length(Phi1)-1);
        R = [ROneSize, Rtmp];
        B1All = 1.01:0.01:10;
        [Bn,RBn,A2,Cp,kk] = FlowSolver(Phi1,R,B, B1All,0);
        C(j) = max(Cp(:));
    end
    EffBlock1(:,k) = 1 - sqrt(16./(27*C));
end
toc

levels = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8];
levelsCols = ["r";"#A2142F";"#D95319";"b";"#0072BD";"c";"g";"#77AC30";"k"];
for jj=1:length(levels)
    hold on
    [C,h] = contour(VelAv,RAlls,EffBlock1,[levels(jj) levels(jj)]);
    
    h.LineColor = levelsCols(jj);
    clabel(C,h,'FontSize',15)
    hold off
end
xlabel('$\phi$','Interpreter','latex')
ylabel('$r$','Interpreter','latex')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Figure4d.eps','Resolution',1200)


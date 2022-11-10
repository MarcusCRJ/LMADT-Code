%% Code for figure 9a

Tol = 0.000001;
B=0.2;
m=0.6;

figure('DefaultAxesFontSize',19,'Position', [10 10 900 600])

NumberStreamtubes = 5;
NumberOfRows = 4;
B1All = 1.01:0.01:10;
KAll = 0.1:0.1:16;
%Fix One Bypass
for i=1:length(KAll)
    Phi1 = ones(NumberStreamtubes,1)';
    R = [1 ; (ones(length(Phi1)-1,1))*(1/B - 1)/(length(Phi1)-1)]';
    Kaim = KAll(i);
    for j=1:NumberOfRows
        [Bn,RBn,A2,Cp,k] = FixedKFinder(Phi1,R,B,Tol,Kaim);
        CpTmp(i,j) = Cp.*((B*sum(R.*Phi1)).^3);
        Uav = B*sum(R.*Phi1);
        if j>1
            BnNew = Bn(2:3)'*RBn(2:3)/(sum(RBn(2:3)));
            Bn = [Bn(1);BnNew;Bn(4:end)];
            RBn = [RBn(1);sum(RBn(2:3));RBn(4:end)];
        end
        R = RBn';
        Phi1=(1-m)*Bn' + m*(Uav);
    end
end
plot(KAll,CpTmp(:,4),'bo-')
xlim([0,16])
hold on
%Fix Two Bypass
for i=1:length(KAll)
    Phi1 = ones(NumberStreamtubes,1)';
    R = [1 ; (ones(length(Phi1)-1,1))*(1/B - 1)/(length(Phi1)-1)]';
    Kaim = KAll(i);
    for j=1:NumberOfRows
        [Bn,RBn,A2,Cp,k] = FixedKFinder(Phi1,R,B,Tol,Kaim);
        CpTmp(i,j) = Cp.*((B*sum(R.*Phi1)).^3);
        Uav = B*sum(R.*Phi1);
        if j>2
            BnNew = Bn(2:3)'*RBn(2:3)/(sum(RBn(2:3)));
            Bn = [Bn(1);BnNew;Bn(4:end)];
            RBn = [RBn(1);sum(RBn(2:3));RBn(4:end)];
        end
        R = RBn';
        Phi1=(1-m)*Bn' + m*(Uav);
    end
end
plot(KAll,CpTmp(:,4),'r-*')
xlim([0,16])
%Fix Three Bypass
for i=1:length(KAll)
    Phi1 = ones(NumberStreamtubes,1)';
    R = [1 ; (ones(length(Phi1)-1,1))*(1/B - 1)/(length(Phi1)-1)]';
    Kaim = KAll(i);
    for j=1:NumberOfRows
        [Bn,RBn,A2,Cp,k] = FixedKFinder(Phi1,R,B,Tol,Kaim);
        CpTmp(i,j) = Cp.*((B*sum(R.*Phi1)).^3);
        Uav = B*sum(R.*Phi1);
        if j>3
            BnNew = Bn(2:3)'*RBn(2:3)/(sum(RBn(2:3)));
            Bn = [Bn(1);BnNew;Bn(4:end)];
            RBn = [RBn(1);sum(RBn(2:3));RBn(4:end)];
        end
        R = RBn';
        Phi1=(1-m)*Bn' + m*(Uav);
    end
end
plot(KAll,CpTmp(:,4),'k')
xlim([0,16])
xlabel('$k$','Interpreter','latex')
ylabel('$C_{P}^{*}$','Interpreter','latex')
legend('One Bypass', 'Two Bypass', 'Three Bypass','Location','southeast')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig9a.eps','Resolution',1200)
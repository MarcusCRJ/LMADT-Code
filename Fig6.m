% Shear Through Turbine Code
% Figures 6a, 6b code.

clear all
tic
TOL = 0.000001;
res = 1000;
B = 0.2;
%B=0.1;
ResThroughTurb=20;

DeltaAll = -0.1:0.05:0.55;
%DeltaAll = -0.05:0.025:0.55;

A2All = 0.03 :0.01:0.98;
%A2All = 0.4 :0.01:0.98;

AllCpShear = zeros(length(A2All),length(DeltaAll));
AllCpUni = AllCpShear;

for j=1:length(DeltaAll)
    delta = DeltaAll(j);
    
    a2av = 1;
    Speed = @(x) ((1-delta) + (x)*2*delta./a2av);
    
    for i=1:length(A2All)
        
        a2av = A2All(i);
        BypassRes = ceil(10/(a2av*B));
        TopSpeed = (1-delta) + (1/B)*2*delta/a2av;
        
        
        RA = a2av*ones(ResThroughTurb,1)/ResThroughTurb;
        RB = (1/B - a2av)*ones(BypassRes,1)/BypassRes;
        
        R = [RA;RB];
        RC = cumsum(R);
        RMid = [0;RC(1:end-1)]./2 + RC./2;
        Phi = Speed(RMid);
        
        [bnall,rnall,CP,k,Ct,UCube,USquare] = FlowSolverVaryNew(R,Phi,B,a2av,TOL,res);
        
        
        CpShear(i) = CP;
        CpShearS(i) = CP/((B*(Phi'*R))^3);
        
        CpShearDraper(i) = CP/UCube;
        CTAllS(i) = Ct/USquare;
        
        % Uniform Through Turbine
        Phi = Speed(RMid);
        
        %% Uniform
        %% Average Inflow
        RTrack = cumsum(R);
        MassTrack = cumsum(R.*Phi);
        ind = find(MassTrack>a2av,1,'first');
        M = ind;
        % Keep a small amount of the bypass to get split by the turbine.
        padding = 0.0001;
        RNew = (a2av - (MassTrack(M-1)))/Phi(M) + padding;
        
        if RNew>R(M)
            R = [R(1:M-1);RNew;R(M+1)+R(M)-RNew;R(M+2:end)];
            Phi = [Phi(1:M-1);(Phi(M)*R(M) + (R(M+1)+R(M)-RNew)*Phi(M+1))/RNew;Phi(M+1:end)];
        else
            Phi = [Phi(1:M-1);Phi(M);Phi(M:end)];
            R = [R(1:M-1);RNew;R(M)-RNew;R(M+1:end)];
            
        end
        
        %% Save Original Profile
        
        Phi(1:M) = Phi(1:M)'*R(1:M) / sum(R(1:M));
        
        
        Phi = [Phi(1);Phi(M+1:end)];
        R = [sum(R(1:M));R(M+1:end)];
        
        
        %% Run uniform for new case
        left = Phi(1,1)*1.000001; right = 25*Phi(2,1); delta = 1e6; %bmplus1=0;
        while delta>TOL
            B1All = [left:(right-left)/res:right];
            [Bn,RBn,A2,Cp,k] = FlowSolver(Phi',R',B,B1All);
            tt=0;
            for ii=2:length(B1All)
                if (isnan(A2(ii-1)) || A2(ii-1)>a2av) && A2(ii)<a2av
                    index = ii;
                    delta = abs(A2(ii)-a2av);
                    left=B1All(index-1); right = B1All(index);
                    tt=1;
                    break
                end
            end
            if tt==0
                [~,midind] = min(abs(A2-a2av));
                delta = abs(left-right);
                left = B1All(midind-1); right = B1All(midind+1);
                
            end
            
        end
        Bn = Bn(:,ii);
        RBn = RBn(:,ii);
        CP = Cp(ii);
        A2 = A2(ii);
        k=k(ii);
        
        
        KAllUni(i) = k;
        CpUni(i) = CP;
        
        
        UStar = Phi(1)^3;
        UStar2 = Phi(1)^2;
        
        CpDraper(i) = CP/UStar;
        CTAllU(i) = k*a2av^2/UStar2;
        CpUniS(i) = CP/((B*(Phi'*R))^3);
        %toc
    end
    AllCpShear(:,j) = CpShearDraper;
    AllCpUni(:,j) = CpDraper;
    
    display(j)
    [MaxShear(j), indshear] = max(CpShearDraper);
    [MaxUni(j), induni] = max(CpDraper);
    
    A2Shear(j) = A2All(indshear);
    A2Uni(j) = A2All(induni);
end

toc

figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
hold on
plot(DeltaAll,1-((16/27)*(1./MaxShear)).^(1/2),'k','LineWidth',2)
plot(DeltaAll,1-((16/27)*(1./MaxUni)).^(1/2),'r','LineWidth',2)
legend('Sheared Flow','Uniform Flow')
xlabel('$\delta$','Interpreter','latex')
ylabel('$B_{\textit{eff}}$','Interpreter','latex')
xlim([-0.05 0.5])
ylim([0 0.7])
%legend('Sheared Flow','Uniform Flow')

exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig6b.eps','Resolution',1200)

hold off
figure('DefaultAxesFontSize',19,'Position', [10 10 700 600])
hold on
plot(DeltaAll,A2Shear,'k')
plot(DeltaAll,A2Uni,'r')
legend('Sheared Flow','Uniform Flow')
xlabel('$\delta$','Interpreter','latex')
ylabel('$B_{eff}$','Interpreter','latex')

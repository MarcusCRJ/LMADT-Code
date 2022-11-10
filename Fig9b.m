%% Code for figure 9b
B = 0.2;
% Switch staggered on or off to get fig 5a or 5b
staggered=1;
rlen=4;
figure('DefaultAxesFontSize',19,'Position', [10 10 900 600])
R = [1.1,1,0.8,1,1.1];
Phi = flip([1,1,1,1,1]);
kall = 0.1:0.25:16;
CpTrack = zeros(length(kall),rlen,length(Bs),length(Ms));
A2All = CpTrack;
lins=["b-o","r-*","k-"];

for tt=2:4
InstantMix = 1;
StreamtubesTotal = tt;
Phi = ones(1,StreamtubesTotal);
rest = (1/B - 1.5)/(StreamtubesTotal-1)*ones(1,StreamtubesTotal-1);
R = [1.5,rest];

BnAll = zeros(StreamtubesTotal+1,length(kall),rlen,length(Bs),length(Ms));
RBnAll = BnAll;

tic
for i=1:length(kall)
    BnTmp = zeros(StreamtubesTotal+1,rlen,length(Bs),length(Ms));
    RBnTmp = zeros(StreamtubesTotal+1,rlen,length(Bs),length(Ms));
    CpTmp = zeros(rlen,length(Bs),length(Ms));
    A2Tmp = zeros(rlen,length(Bs),length(Ms));
    for ii=1:length(Bs)
        B=Bs(ii);
        for jj=1:length(Ms)
            m=Ms(jj);
            PhiIn = Phi';
            R = (1/B - 1.5)*ones(1,length(Phi)-1)/(length(Phi)-1);
            RIn = [1.5, R]';
            Cp=0;
            for j=1:rlen
                if isnan(Cp)
                CpTmp(j,ii,jj) = Cp;
                A2Tmp(j,ii,jj) = NaN;
                %BnAll(1:length(Bn),i,j,ii,jj) = Bn;
                BnTmp(1:length(Bn),j,ii,jj) = NaN;
                %RBnAll(1:length(Bn),i,j,ii,jj) = RBn;
                RBnTmp(1:length(RBn),j,ii,jj) = NaN;                      
                else
                    
                [Bn,RBn,A2,Cp,~] = FixedKFinder(PhiIn',RIn',B,0.000001,kall(i));
                CpTmp(j,ii,jj) = Cp;
                A2Tmp(j,ii,jj) = A2;
                BnTmp(1:length(Bn),j,ii,jj) = Bn;
                RBnTmp(1:length(RBn),j,ii,jj) = RBn;
                if InstantMix == 1
                    bb=2;
                    while bb<length(Bn) && length(Bn)>StreamtubesTotal
                       if Bn(bb)==Bn(bb+1)
                           RBn(bb+1) = RBn(bb)+RBn(bb+1);
                           Bn = [Bn(1:bb-1);Bn(bb+1:end)];
                           RBn = [RBn(1:bb-1);RBn(bb+1:end)];
                           bb=bb-1;
                       end
                       bb = bb+1;
                    end
                    if StreamtubesTotal == 2 && length(Bn)>StreamtubesTotal
                       Bprime = (Bn(2)*RBn(2)+Bn(3)*RBn(3))/(RBn(2)+RBn(3));
                       Bn = [Bn(1),Bprime]';
                       RBn = [RBn(1),RBn(2)+RBn(3)]' ;
                    else
                    while length(Bn)>StreamtubesTotal
                        Bprime = (Bn(2)*RBn(2)+Bn(3)*RBn(3))/(RBn(2)+RBn(3));
                        Bn = [Bn(1),Bprime,Bn(4:end)']';
                        RBn = [RBn(1),RBn(2)+RBn(3),RBn(4:end)']';
                    end
                    end
                end
                PhiIn = m*(sum(Bn.*RBn)*B) + (1-m)*Bn;
                
                RIn = RBn;
                if staggered ==1
                    PhiIn(RIn==0)=0;
                    RIn(PhiIn==0)=0;
                    PhiTmp = PhiIn(PhiIn>0);
                    PhiIn = zeros(size(PhiTmp));
                    PhiIn(1:sum(PhiTmp>0)) = flip(PhiTmp(PhiTmp>0));
                    RTmp = RIn;
                    RIn = zeros(size(RTmp));
                    RIn(1:sum(RTmp>0)) = flip(RTmp(RTmp>0));
                end
                PhiIn = PhiIn./(Phi(1));
                end
            end
        end
    end
    BnAll(:,i,:,:,:)=BnTmp;
    RBnAll(:,i,:,:,:)=RBnTmp;   
    CpTrack(i,:,:,:) = CpTmp;
    A2All(i,:,:,:) = A2Tmp;
end
toc
TurbineNumber = 4;
plot(kall,CpTrack(:,TurbineNumber,1,1),lins(tt-1))
hold on

end
xlabel('$k$','Interpreter','latex')
ylabel('$C_{P}^{*}$','Interpreter','latex')
legend('One Bypass', 'Two Bypass', 'Three Bypass','Location','southeast')
exportgraphics(gcf,'C:/Users/Marcus/Downloads/Fig9b.eps','Resolution',1200)
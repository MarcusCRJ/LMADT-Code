%% Code for figure 12a 
%Resolution of operating conditions checked.
Res = 150; %How many B1 points to check.


%Putting on matlabs parallel option, need to replace for with parfor for
%first loop.
%parpool('local',16)

%List of blockages and mixing rates.
Bs = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13, ...
    0.14,0.15,0.16,0.17,0.18,0.19,0.20];
Ms = [0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72, ...
    0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85, ... 
    0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99];

simSpace = [length(Bs),length(Ms)];
numSims = prod(simSpace);
GlobalMaximums = zeros(numSims,1);
k1 = GlobalMaximums;
k2=k1;
k3=k1;
tic

for idx=1:numSims
% Set Inputs
scaletop = 1.75;
Bs = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.22,0.24,0.26,0.28,0.3];
Ms = [0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99];

[ind1, ind2] = ind2sub(simSpace, idx);
B1 = logspace(log10(1.0001),log10(scaletop),Res);
B=Bs(ind1);
Phi = [1,1,1,1,1];
Rtmp = (1/B - 1.5)*ones(1,length(Phi)-1)/(length(Phi)-1);
R = [1.5, Rtmp];
m=Ms(ind2); %Mixing Rate
BnVals = length(Phi)+4;
%(BnValue, Turb1_B1Idx, Turb2_B1Idx, Turb3_B1Idx)
BnAll = zeros(BnVals,Res+1,Res+1,Res+1);
PhiAll = BnAll;
RBnAll = BnAll;
CpAll = zeros(Res+1,Res+1,Res+1); %(Turb1_B1Idx, Turb2_B1Idx, Turb3_B1Idx)
CpTrack = zeros(Res+1,Res+1,Res+1);
KAll = CpTrack;

%No turbines at all
BnAll(1:length(Phi),1,1,1)=Phi;
PhiAll(1:length(Phi),1,1,1)=Phi;
RBnAll(1:length(Phi),1,1,1)=R;

%Flow through only first turbine all others non-operational
[Bn,RBn,A2,Cp,k,~,~] = FlowSolver(Phi,R,B,B1);

BnAll(1:length(Bn(:,1)),2:end,1,1)=Bn;%./Bn(1,:);
PhiAll(1:length(Bn(:,1)),2:end,1,1)=(1-m)*Bn +(m)*(sum(Bn.*RBn,1)*B);
RBnAll(1:length(Bn(:,1)),2:end,1,1)=RBn;
CpAll(2:end,1,1)= Cp;
CpTrack(2:end,1,1) = Cp./3;
KAll(2:end,1,1)=k;

%Loop through first and second turbine
for ii=2:Res+1
    if isnan(CpAll(ii,1,1))
        BnAll(1:length(Bn(:,1)),ii,2:end,1)=NaN;
        PhiAll(1:length(Bn(:,1)),ii,2:end,1)=NaN;
        RBnAll(1:length(Bn(:,1)),ii,2:end,1)=NaN;
        CpAll(ii,2:end,1)= NaN;
        CpTrack(ii,2:end,1)= NaN;
        KAll(ii,2:end,1)=NaN;
    else
        PhiTmp = PhiAll(1:length(Phi)+1,ii,1,1)';
        PhiIn = PhiTmp;
        RTmp = RBnAll(1:length(Phi)+1,ii,1,1)';
        RIn = RTmp;
        [Bn,RBn,A2,Cp,k,~,~] = FlowSolver(PhiIn,RIn,B, logspace(log10(1.0001*PhiIn(1)),log10(scaletop*PhiIn(1)),Res));
        BnAll(1:length(Bn(:,1)),ii,2:end,1)=Bn;
        PhiAll(1:length(Bn(:,1)),ii,2:end,1)=(1-m)*Bn +(m)*(sum(Bn.*RBn,1)*B);
        RBnAll(1:length(Bn(:,1)),ii,2:end,1)=RBn;
        CpAll(ii,2:end,1)= Cp;
        CpTrack(ii,2:end,1) = (Cp./3) + CpTrack(ii,1,1);
        KAll(ii,2:end,1)=k;
    end
end
%Loop through first and second and third turbine
for ii=2:Res+1
    for jj=2:Res+1
        if isnan(CpAll(ii,jj,1))
            BnAll(1:length(Bn(:,1)),ii,jj,2:end)=NaN;
            PhiAll(1:length(Bn(:,1)),ii,jj,2:end)=NaN;
            RBnAll(1:length(Bn(:,1)),ii,jj,2:end)=NaN;
            CpAll(ii,jj,2:end)= NaN;
            CpTrack(ii,jj,2:end)=NaN;
            KAll(ii,jj,2:end)=NaN;
        else
            PhiTmp = PhiAll(1:length(Phi)+2,ii,jj,1)';
            PhiIn = PhiTmp;
            
            RTmp = RBnAll(1:length(Phi)+2,ii,jj,1)';
            RIn = RTmp;
            [Bn,RBn,A2,Cp,k,~,~] = FlowSolver(PhiIn,RIn,B, logspace(log10(1.0001*PhiIn(1)),log10(scaletop*PhiIn(1)),Res));
            BnAll(1:length(Bn(:,1)),ii,jj,2:end)=Bn;
            PhiAll(1:length(Bn(:,1)),ii,jj,2:end)=(1-m)*Bn +(m)*(sum(Bn.*RBn,1)*B);
            RBnAll(1:length(Bn(:,1)),ii,jj,2:end)=RBn;
            CpTrack(ii,jj,2:end) = Cp./3 + CpTrack(ii,jj,1);
            KAll(ii,jj,2:end)=k;
            
        end
    end
end

[GlobalMax,b]=max(CpTrack(:));
[l,m,n,o]=ind2sub(size(CpTrack),b);
GlobalMaximums(idx) = GlobalMax;
k1(idx)=KAll(l,1,1);
k2(idx)=KAll(l,m,1);
k3(idx)=KAll(l,m,n);
end 
GlobalMaximums = reshape(GlobalMaximums, simSpace);
k1 = reshape(k1, simSpace);    
k2 = reshape(k2, simSpace);
k3 = reshape(k3, simSpace);
toc

save('align_variable_k_3turb','GlobalMaximums','k1','k2','k3','k4','Bs','Ms')

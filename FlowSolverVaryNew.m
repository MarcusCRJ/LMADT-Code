function [Bn,RBn,Cp,k,Ct,UCube,USquare] = FlowSolverVaryNew(R,Phi,B,a2av,TOL,res)
%% Average Inflow
RTrack = cumsum(R);
MassTrack = cumsum(R.*Phi);
ind = find(MassTrack>a2av,1,'first');
M = ind;
while M==1
    R = [R(1)/2;R(1)/2;R(2:end)];
    Phi = [Phi(1);Phi(1);Phi(2:end)];
    RTrack = cumsum(R);
    MassTrack = cumsum(R.*Phi);
    ind = find(MassTrack>a2av,1,'first');
    M = ind;
end
% Keep a small amount of the bypass to get split by the turbine.
padding = 0.00001;

RNew = (a2av - (MassTrack(M-1)))/Phi(M) + padding;


if RNew>R(M)
    R = [R(1:M-1);RNew;R(M+1)+R(M)-RNew;R(M+2:end)];
    Phi = [Phi(1:M-1);(Phi(M)*R(M) + (R(M+1)+R(M)-RNew)*Phi(M+1))/RNew;Phi(M+1:end)];
else
   Phi = [Phi(1:M-1);Phi(M);Phi(M:end)];
R = [R(1:M-1);RNew;R(M)-RNew;R(M+1:end)]; 
    
end

USquare = sum(R(1:M).*Phi(1:M).^3)/ sum(R(1:M));
UCube = sum(R(1:M).*Phi(1:M).^2)/ sum(R(1:M));
%% Save Original Profile
PhiFull = Phi;
RFull = R;
Phi(1:M) = Phi(1:M)'*R(1:M) / sum(R(1:M)); 
Phi = [Phi(1);Phi(M+1:end)];
R = [sum(R(1:M));R(M+1:end)];


momcorr = sum(RFull(1:M).*PhiFull(1:M).^2) - R(1)*Phi(1)^2;

%% Run uniform for new case
left = Phi(1,1)*1.000001; right = 25*Phi(1,1); delta = 1e6; %bmplus1=0;
while delta>TOL
    B1All = [left:(right-left)/res:right];
    [Bn,RBn,A2,Cp,k] = FlowSolver(Phi',R',B,B1All,momcorr);
    tt=0;
    for i=2:length(B1All)
        %[Bn,RBn,A2,Cp,k,PhiAv,RAv] = FlowSolver(Phi',R',B,B1All(i),momcorr);
        if (isnan(A2(i-1)) || A2(i-1)>a2av) && A2(i)<a2av
            index = i;
%             if i==1
%                 %Acts as a propellor
%                 delta=0;
%                 Cp=NaN;
%                 break
%             end
            delta = abs(A2(i)-a2av);
            left=B1All(index-1); right = B1All(index);
            tt=1;
            break
        end
    end
    if tt==0
       [~,midind] = min(abs(A2-a2av));
       delta = abs(left-right);
       if midind ==1
          left = B1All(1); right = B1All(3); 
       else
       left = B1All(midind-1); right = B1All(midind+1);
       end
    end

end
Bn = Bn(:,i);
RBn = RBn(:,i);
Cp = Cp(i);
A2 = A2(i);


%% Find Wake
if isnan(A2)
    Cp = NaN;
else
    
left = 0.000001; right = PhiFull(M)*1.1; delta = 1e6; trackval = 0;%bmplus1=0;
while delta>TOL
    B1All = [left:(right-left)/res:right];
    BAll = (B1All.^2 - PhiFull(M)^2 + PhiFull(1:M).^2).^(1/2);
    RBAll = RFull(1:M).*PhiFull(1:M)./BAll;
    tt = 0;
    for i=2:length(B1All)
        Rtrack = sum(RBAll,1);
        if isreal(Rtrack(:,i-1)) && ((Rtrack(:,i-1)<RBn(1) && Rtrack(:,i)>RBn(1))...
                || (Rtrack(:,i-1)>RBn(1) && Rtrack(:,i)<RBn(1)))
            index = i;
            %delta = abs(RBAll(:,i)'*BAll(:,i)-a2av);
            left=B1All(index-1); right = B1All(index);
            delta = abs(right-trackval); trackval = right;
            tt=1;
            break
        end
    end
    
    if tt==0
       tmp = find(imag(Rtrack) ==0,1,'first');
       left = B1All(tmp-1); right = B1All(tmp); 
    end
    
end
end
%% Adjust CP, Obtain K, CT OTHERS???
T = (Cp/a2av + RBn(1)*Bn(1)^2 - sum(RBAll(:,index).*BAll(:,index).^2));
Cp = T*a2av;

Bn = [BAll(:,index);Bn(2:end)];
RBn = [RBAll(:,index);RBn(2:end)];

% Need to guess k:

k = T/a2av^2;
Ct = T;


end
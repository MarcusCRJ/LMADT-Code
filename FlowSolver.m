function [Bn,RBn,Aav,Cp,k] = FlowSolver(Phi,R,B,B1,momcorr)
%% Flow solver for uniform flow through the turbine
% Function takes in:
% 1) vectors Phi and R which are the speed and width of each streamtube.
% 2) scalar B which is the blockage ratio.
% 3) scalar or vector B1 which is \beta_{1} (or range of \beta_{1}'s), the 
% downstream flow speed of the streamtube which most immediately bypasses 
% the turbine.
% 4) optional scalar momcorr which is the momentum correction term used in
% calculating the non-uniform flow through the turbine case.

% Outputs:
% 1) Arrays Bn and RBn which are the downstream streamtube flowspeeds and
% sizes for each chosen B1.
% 2) Vectors Aav, Cp, k are the average flowspeed throught the turbine 
% \alpha_av, the power coefficient C_{P} and resistance coefficient,
% respectively.

% Author: Marcus C. R. Juniper <marcus.juniper@eng.ox.ac.uk>
% Paper: TBC

% If no momentum correction defined set it to 0.
if ~exist('momcorr','var')
    momcorr = 0;
end

% Pre-Allocate array sizes.
Bn=zeros(length(Phi)+1,length(B1));
RBn = Bn;
Cp = zeros(1,length(B1));
Aav = Cp;
k = Cp;


% If any of the flow has null values we set outputs to NaN and skip the
% computation.
if isnan(Phi)
    Bn(:,:) = NaN;
    RBn(:,:) = NaN;
    Cp(:) = NaN;
    Aav(:) = NaN;
    k(:) = NaN;
else
    for i=1:length(Phi)
        Bn(i+1,:) = (B1.^2 - Phi(1)^2 + Phi(i)^2).^(1/2);
    end
    RBn = R(2:end).*Phi(2:end)./Bn(3:end,:)';
    
    C1 = 2*(R(1)*Phi(1)*(Phi(1)-B1') + sum( R(2:end).*Phi(2:end).*(Phi(2:end) - Bn(3:end,:)'),2)) ...
        + (1/B)*(B1'.^2 -Phi(1)^2) - B1'.^2 + 2*momcorr;
    
    C2 = R(1)*Phi(1) - B1'.*(1/B - sum(RBn,2));
    
    % Find roots of the cubic equation.
    B0 = CardanRoots([ones(length(B1),1), -(2*C2 + B1'), (2*C2.*B1' +C1), -C1.*B1'],0)';

    % Find the physically allowable solution.
    B0True = B0(B0<Phi(1) & B0>0);
    B0All = NaN(1,length(B1));
    B0All(logical(sum((B0<Phi(1) & B0>0),1))) = B0True;
    Bn(1,:) = B0All;
    
    Aav = (B0All'.*(C2))./(B0All'-B1');
    
    Bla=(Aav<0);
    Aav(Bla) = NaN;
    
    RBn = [(Aav'./B0All)' ,(1/B - sum([Aav./B0All',RBn],2)), RBn]';
    
    Remove = logical(sum(RBn<0,1));
    RBn(:,Remove)=NaN;
    Bn(:,Remove)=NaN;
    Aav(Remove) = NaN;
    %Cp = ((B1'.^2 - B0All'.^2).*Aav)./(((B*sum(R.*Phi)).^3)'); %CP*
    %calculation
    Cp = ((B1'.^2 - B0All'.^2).*Aav);
    k = (B1'.^2 - B0All'.^2)./(Aav.^2);
end
end
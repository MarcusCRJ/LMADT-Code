function [Bn,RBn,A2,Cp,k] = FixedKFinder(Phi,R,B,Tol,Kaim)


Diff = 1e6; left = Phi(1)*1.0001; right = Phi(1)*2; tmp=0;
Res = 500;
while Diff>Tol
B1 = linspace(left,right,Res);
[Bn,RBn,A2,Cp,k] = FlowSolver(Phi,R,B, B1);

K = k-Kaim;
t=0;
for i=2:length(k)
    if (K(i-1)<0 && K(i)>0) || (K(i-1)>0 && K(i)<0)
        KInd = k(i-1); left = B1(i-1); right = B1(i);
        t=1;
        break
    end
end

if t==0
 right = B1(end-50);
end

Diff = abs(right-tmp); tmp = right;

end


Bn = Bn(:,1);
RBn = RBn(:,1);
A2 = A2(1);
Cp = Cp(1);
k = k(1);
end
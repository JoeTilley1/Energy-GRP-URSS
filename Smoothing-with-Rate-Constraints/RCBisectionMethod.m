function newZ = RCBisectionMethod(Zmax,Zmin,n,tau,M,s)
%RCBISECTIONMETHOD: Given our ideal lambda, we have Zmin and Zmax. If
%one of these is feasible, we keep it and store the times of exit.
%Otherwise, we use bisection method on Zmin and Zmax (which doesnt change
%total variation) to arrive at a level Z which is feasible.

%% Testing if Zmax or Zmin are feasible

Type = RCZtypefun(Zmax,Zmin,n,tau,M,s);

if Type == 4
    Z=Zmax;
elseif Type == 5
    Z=Zmin;
else
    %% Bisection method
    Z=0.5*(Zmax+Zmin);
    [m,~,~]=RCZtypefun2(Z,n,tau,M,s);
    while m==1 || m==3
        if m==1 %(Exits bottom)
            Zmin = Z;
        else %m=3 (Exits top)
            Zmax = Z;
        end
        Z=0.5*(Zmax+Zmin);
        [m,~,~]=RCZtypefun2(Z,n,tau,M,s);
    end
end

%Cutting off when contraint is broken
newZ=[];
for i=1:length(Z)
    newZ(i) = Z(i);
    if newZ(i)<-s || newZ(i)>M+s;
        break
    end
end

end


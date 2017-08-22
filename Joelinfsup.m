function Lam=Joelinfsup(p,tau,X0,L,M,a1,a2,s)
%This file calculates the multiplier L by using a bisection method

Z=JoeLevelfun(p,tau,X0,L,M,a1,a2,s);     % Initial candidate strategy of storage levels
n=numel(p);
format long

%% Adjusting L

%Initial estimates
lmax=1e10;
lmin=-1e10;
l=L;

%Type of constraint that is broken
m=JoeZtypefun(Z,n,tau,M,s);
m=m(1);     

while abs(lmax-lmin)>1e-14 
%% Adjustment of L by bisection
    switch m;
        case {1}
            lmin=l;         %bigger
        case {3}
            lmax=l;         %smaller
        case {2,4,5}
            lmax=l;
            lmin=l;         %Feasible
    end
    l=0.5*(lmin+lmax);
    Z=JoeLevelfun(p,tau,X0,l,M,a1,a2,s);
    m=JoeZtypefun(Z,n,tau,M,s);
    m=m(1);
end

Lam=l;



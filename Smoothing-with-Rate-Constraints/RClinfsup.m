function [Zmin,Zmax,Lam] = RClinfsup(p,tau,X0,M,a1,a2,s,q1,q2)
%This file calculates the multiplier L by increasing or decreasing as appropriate.
%Firstly calculated Ltop and Lbot, the maximum and minimum values of lambda which give feasible strategies.
%Then calculates Lam \in [Ltop,Lbot] such that Lam gives smallest Total Variation.

n=numel(p);
format long

%% Special cases

%Check if the max power is too small
maximum = alpha(q1,a1,a2).*q1-p(tau);
[maxZmax,maxZmin] = RCLevelfun(p,tau,X0,maximum,M,a1,a2,s,q1,q2);
maxm = RCZtypefun(maxZmax,maxZmin,n,tau,M,s);
if maxm == 1
    disp('We have a special case where we need max power and thats too small. Break the programme.')
end

%Check if the min power is too big
minimum = alpha(-q2,a1,a2).*-q2-p(tau);
[minZmax,minZmin] = RCLevelfun(p,tau,X0,minimum,M,a1,a2,s,q1,q2);
minm = RCZtypefun(minZmax,minZmin,n,tau,M,s);
if minm == 3
    disp('We have a special case where we need min power and thats too much. Break the programme.')
end



%% Determining Ltop (the largest value of L which gives a feasible strategy)
topmax= alpha(q1,a1,a2).*q1-p(tau);
topmin= alpha(-q2,a1,a2).*-q2-p(tau);
feas=0; %to enter while loop

while ~((abs(topmax-topmin)<1e-10) && (feas==2 || feas==4 || feas==5)) %Make Ltop precise enough and also feasible
    % Adjustment of Ltop
    [Zmax,Zmin] = RCLevelfun(p,tau,X0,0.5*(topmax+topmin),M,a1,a2,s,q1,q2);
    m = RCZtypefun(Zmax,Zmin,n,tau,M,s);
    switch m;
        case {1,2,4,5}                          %average is feasible or too small
            topmin=0.5*(topmax+topmin);         
            feas = m;                           %break loop only if topmin is feasible
        case {3}                                %average is too large
            topmax=0.5*(topmax+topmin);         
    end
end

Ltop=topmin;

%% Determining Lbot (the smallest value of L which gives a feasible strategy)
botmax= alpha(q1,a1,a2).*q1-p(tau);
botmin= alpha(-q2,a1,a2).*-q2-p(tau);
feas=0; %Initial value to enter while loop

while ~( (abs(botmax-botmin)<1e-10) && (feas==2 || feas==4 || feas==5) ) %Make Lbot precise enough and also feasible
    % Adjustment of Lbot
    [Zmax,Zmin] = RCLevelfun(p,tau,X0,0.5*(botmax+botmin),M,a1,a2,s,q1,q2);
    m = RCZtypefun(Zmax,Zmin,n,tau,M,s);
    switch m;
        case {2,3,4,5}                          %average is feasible or too large
            botmax=0.5*(botmax+botmin);         
            feas = m;                           %break loop only if botmax is feasible
        case {1}                                %average is too small
            botmin=0.5*(botmax+botmin);         
    end
end

Lbot=botmax;


%% We have found the max and min lambdas in our range of feasibility
% We now calculate the total variation for each lambda in this range and
% choose the smallest.
% There is another method (which may be a computationally
% quicker) by determining the L in Lrange which takes the longest time before 
% hitting rate constraint and making this Lam.

Lrange = linspace(Lbot,Ltop,ceil((Ltop-Lbot)/s)); %Determines Lambda to within s

for i=1:length(Lrange)
    x=[];           %list storing values of x
    y=[];           %list storing values of y
    Lvec(1) = Lrange(i);    %lambda_1
    % Complete the method of finding the net energy over time (as in RCLevelfun.m)
    for k=1:n+1-tau
        alphaxx = p(k+tau-1)+Lvec(k);   %ideal value of alpha(x(k))*x(k) but may be x(k) \notin X
        if alphaxx <= -a2*q2
            x(k) = -q2;
        elseif alphaxx >= q1/a1
            x(k) = q1;
        elseif p(k+tau-1)+Lvec(k)>=0
            x(k) = a1*(p(k+tau-1)+Lvec(k));
        else
            x(k) = (p(k+tau-1)+Lvec(k))/a2;
        end
        y(k) = alpha(x(k),a1,a2)*x(k)-p(k+tau-1);
        Lvec(k+1) = y(k);
    end
    % Calculate Total Variation
    TV(i)=0;
    for k=1:length(y)-1
        TV(i)=TV(i)+abs(y(k+1)-y(k));
    end
end

%Finding the minimum Total Variation and corresponding Lambda and Levels
[~,idealindex] = min(TV);
Lam = Lrange(idealindex);
[Zmax,Zmin] = RCLevelfun(p,tau,X0,Lam,M,a1,a2,s,q1,q2);

end


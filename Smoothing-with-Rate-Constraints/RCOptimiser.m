function y=RCOptimiser(M,a1,a2,Q1,Q2,s)

%% Objective
%To minimise the total variation of the net wind energy, assuming that the store has no rate constraints.
%The (negative of the) net wind energy at time t is f(p_t,x_t)=alpha(x_t)*x_t-w_t, where w_t is wind energy at time t , x_t is the energy taken off the grid (or put onto if negative) if the store charges up (or discharges if -ve) with x_t units of energy, and alpha(.) the ineffiency.  
%The objective is equivalent to stretching a string i.e. trying to make f as straight as possible.  The gradient of this string satisfies the same properties as the multiplier in our usual price arbitrage problem.
%This means that we should look for a multipler (L_1,...,L_n) in the usual way and then set x_t satisfying x_t=(L_t+w_t)/(alpha(x_t)) for all times t.

%One timestep in this program is defined by the wind energy series - if the series is half-hourly, then times t=1,2,...,n correspond to half hours.  
%We solve min_x TV(sum_t a(x_t)*x_t-p_t), where TV=total variation.

%Wind energy can be replaced with any other time series we wish to smooth using a store e.g. supply-demand.

%% Assumptions and notation
%We assume that the store starts at level X0 at time 0 and finishes empty at time n (X0 currently set at 0).
%At time t>=1, the store's cost is related to the price p_t.  The store acts accordingly over the interval (t-1,t].
%Time t is called an "update time" if the multiplier changes value at that time.  In this program, we calculate times tau - this always refers to an update time + 1 timestep.
%Time T is called an "exit time" if we need to look no further than the price p_T in order to decide how to update the multiplier.

%% Inputs
%Need to input a wind power series p.  (Done here via a csv file.) 

%Input "time" - this is the timelength in hours of one timestep in the wind power series. 

%Energy capacity units should be consistent with the units of the data series

%M>0 is energy capacity of the store

%The following are dimensionless:

%a1 \in (0,1] is charge efficiency
%a2 \in (0,1] is discharge efficiency

%s>0 is allowed error i.e. we say that the store is full if the level of stored energy lies in (M-s,M] and it is empty if the level lies in [0,s).

%Q1 and Q2 are the rate constraints.
%Their units are consistent with M (e.g. M in MWh, Q1 and Q2 in MW).
%Then x_t \in [-Q2*time, Q1*time] where "time" is one timestep.

%% Outputs

%"Level": The optimal storage level at each time.  The units are consistent with the input energy capacity and rates (M,Q1,Q2). 
%"Rates": the charge and discharge powers required to achieve the maximum level of smoothing given the capacity constraints.
%"TotalVariation": the minumum achievable total variation of the store's output minus wind energy (i.e. the net energy)

tic

%% Input wind series p

%p=csvread('Irishwind.csv');            %About a month's worth of quarterly-hourly wind power in MW.
%p=csvread('Elexonwind2014.csv');       %A year's worth of half-hourly wind power 2014 in MW.
p=csvread('19Jul16AnnualPVData.csv');   %A year's worth of half-hourly university solar power from 19Jul16-19Jul17 in kWh
time=1/2;                             %This is the length of one timestep (in hours) in the wind power series.  ADJUST THIS AS NEEDED. 
%p=p*time;                             %Converting power into energy per timestep (FOR IRISH AND ELEXON ONLY) 
p=p*10^(-3);                          %Converting from kWh to MWh (FOR UNIVERSITY DATA ONLY)
p=transpose(p);
n=numel(p);

%% Converting max power to max power per timestep
q1=Q1*time;
q2=Q2*time;

%% Initial inputs
tau=1;    %Update time +1
X0=0;     %Start level of the store - can change this


%% Keeping track of the variables
Lvec=[];        %List of multipliers (values of lambda_1 used)
tauvec=tau;     %List of (update times+1), starting at time 1   
Level=X0;       %Level of stored energy at the end of each time, starting at time 0                
Tvec=[];        %List of (exit times+1)
X0vec=X0;       %List of energy levels at update times  


while tau<n+1 && -s<=X0 && X0<=M+s;

%% Finding the multiplier and the strategy from time tau to the exit time.  Here tau is the update time +1.

[Zmin,Zmax,L]=RClinfsup(p,tau,X0,M,a1,a2,s,q1,q2);  %Actual multiplier and max/min levels
Z = RCBisectionMethod(Zmax,Zmin,n,tau,M,s); %Optimised store levels until exit time 

%The last entry of Z is the first time that a capacity constraint would be broken if we did not update the multiplier L before then

sizez=numel(Z);  

[Ztype,ttop,tbottom]=RCZtypefun2(Z,n,tau,M,s); %Characterises the shape of the graph of Z

%% New update time
T=tau+sizez;     %Exit time +1
tau1=tau;

if Ztype==2;
    tau=tau1+ttop;        %The time at which the store last hits the top before exiting from the bottom
elseif Ztype==4;
    tau=tau1+tbottom;        %The time at which the store last hits the bottom before exiting from the top
else
    tau=tau1+sizez;
end
    
%The store is either full or empty at the update time (tau-1). 

%% Level of store at the update time (tau-1)
X0=Z(tau-tau1);

%% Capacity checks: finish early if capacity constraints are broken

Capacitycheck=0;
if X0<-s || X0>M+s;
    Capacitycheck=1
    break
end

%% Truncated strategy until new update time (tau-1)
Zc=Z*eye(sizez,tau-tau1);  %Zc includes the level at the update time (tau-1) as the final entry. 
Level=[Level,Zc];    %The ith element of Level is the level at time i-1 (with the first entry being time 0).


%% Keeping track of the running variables
Lvec=[Lvec,L];                  %Multipliers
tauvec=[tauvec,tau];            %Update times +1
Tvec=[Tvec,T];                  %Exit times +1
X0vec=[X0vec,X0];               %Storage levels at update times

end

%We have now found the optimal strategy, with Level=level of store at each time (including time 0).

%% Plots 
x=zeros(1,n);
for i=1:n;
    if Level(i+1)-Level(i)>0;
        x(i)=(Level(i+1)-Level(i))/a1;
    else
        x(i)=(Level(i+1)-Level(i))*a2;
    end
end


%Level of the store
z=linspace(1,n+1,n+1);
figure
plot((z-1)*time,Level(z),'k','Linewidth',2);   
title('Level of the store over time period');
xlabel('Time (h)')
ylabel('Store level (kWh)') %Take care to check units are correct
axis([0 inf -s M+s]);

%Wind energy
z=linspace(1,n,n);
figure;
plot(z*time,p(z),'k','Linewidth',2);
title('Solar energy output over time period');
xlabel('Time (h)')
ylabel('Solar energy (MWh)')
axis([0 inf 0 0.06]);

%Net energy
figure
plot(z*time,-(x(z)-p(z)),'k','Linewidth',2);      
title('Net energy output over time period');
xlabel('Time (h)')
ylabel('Net energy (kWh)')
axis([0 inf 0 max(p)]);


%% Outputs
 Charge=max(x)/(time);                    
 Discharge=min(x)/(time);

 Rates=[Charge,Discharge]   %This is the minimum charge and discharge powers required to smooth wind as much as the capacity of the store allows.

 
 xr=zeros(1,n-1);
 pr=xr;
 
 for i=1:n-1;
     xr(i)=x(i+1);
     pr(i)=p(i+1);
 end
 

z=linspace(1,n-1,n-1); 
TotalVariation=sum(abs(xr(z)-pr(z)-x(z)+p(z)))   %The required total variation corresponding to our optimal strategy.
TotalVariationwind=sum(abs(-pr(z)+p(z)))        %The total variation due of the wind alone

toc

end


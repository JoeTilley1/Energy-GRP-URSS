function Ztypevec=JoeZtypefun(Z,n,tau,M,s)
% Given a level over time Z, we calculate the type of contraint that is broken and times

sizez=numel(Z);

%% Keeping track of when the store hits the top or bottom

Z2=zeros(1,sizez-1);   %Used to keep track if the store hits the top or the bottom before the exit time

for i=1:sizez-1;
        if abs(Z(i)-M)<=s;   %Hits top
            Z2(i)=1;
        elseif abs(Z(i))<=s;   %Hits bottom
            Z2(i)=-1;
        else
            Z2(i)=0;
        end
end

%ttop is the last time Z hits the top; tbottom the last time it hits the bottom (excluding the final time of Z)
ttop=sizez;
tbottom=ttop;

if max(Z2)==1;
    ttop=find(Z2==1,1,'last');
else
    ttop=sizez;
end

if min(Z2)==-1;
    tbottom=find(Z2==-1,1,'last');
else
    tbottom=sizez;
end

%% Defining the types

if Z(sizez)<-s && max(Z)<M-s;       %Type 1: Exits bottom and doesn't hit top beforehand (increase L)
    Ztype=1;
elseif Z(sizez)<=s && max(Z)>=M-s;  %Type 2: Hits or exits bottom, and hits top beforehand (keep L)
    Ztype=2;
elseif Z(sizez)>s && min(Z)>s;      %Type 3: Exits top and doesn't hit bottom beforehand (decrease L)
    Ztype=3;
elseif Z(sizez)>s && min(Z)<=s      %Type 4: Exits top and hits bottom beforehand (keep L)
    Ztype=4;
elseif sizez==n+1-tau && tbottom<sizez; %Type 4: Feasible but hits the bottom before the end time (keep L)
    Ztype=4;
else
    Ztype=5;                       %Either feasible or will need adjusting
end

Ztypevec=[Ztype,ttop,tbottom];


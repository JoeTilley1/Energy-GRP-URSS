function Ztype = RCZtypefun(Zmax,Zmin,n,tau,M,s)
%Function calculates the type of constraints that are broken by the levels
%Zmax and Zmin

%Cutting off when contraint is broken for max
for i=1:length(Zmax)
    newZmax(i) = Zmax(i);
    if newZmax(i)<-s || newZmax(i)>M+s;
        break
    end
end

%Cutting off when contraint is broken for min
for i=1:length(Zmin)
    newZmin(i) = Zmin(i);
    if newZmin(i)<-s || newZmin(i)>M+s;
        break
    end
end

sizezmax=numel(newZmax);
Zmaxtype=0;
sizezmin=numel(newZmin);
Zmintype=0;

%% Keeping track of when the max store hits the top or bottom
Zmax2=zeros(1,sizezmax-1);   %Used to keep track if the store hits the top or the bottom before the exit time

for i=1:sizezmax-1;
    if abs(newZmax(i)-M)<=s;   %Hits top
        Zmax2(i)=1;
    elseif abs(newZmax(i))<=s;   %Hits bottom
        Zmax2(i)=-1;
    else
        Zmax2(i)=0;
    end
end

%ttopmax is the last time newZmax hits the top; tbottommax the last time it hits the bottom (excluding the final time)
ttopmax=sizezmax;
tbottommax=ttopmax;

if max(Zmax2)==1;
    ttopmax=find(Zmax2==1,1,'last');
else
    ttopmax=sizezmax;
end

if min(Zmax2)==-1;
    tbottommax=find(Zmax2==-1,1,'last');
else
    tbottommax=sizezmax;
end


%% Keeping track of when the min store hits the top or bottom
Zmin2=zeros(1,sizezmin-1);   %Used to keep track if the store hits the top or the bottom before the exit time

for i=1:sizezmin-1;
        if abs(newZmin(i)-M)<=s;   %Hits top
            Zmin2(i)=1;
        elseif abs(newZmin(i))<=s;   %Hits bottom
            Zmin2(i)=-1;
        else
            Zmin2(i)=0;
        end
end

%ttop is the last time Z hits the top; tbottom the last time it hits the bottom (excluding the final time of Z)
ttopmin=sizezmin;
tbottommin=ttopmin;

if max(Zmin2)==1;
    ttopmin=find(Zmin2==1,1,'last');
else
    ttopmin=sizezmin;
end

if min(Zmin2)==-1;
    tbottommin=find(Zmin2==-1,1,'last');
else
    tbottommin=sizezmin;
end


%% Defining the types

if newZmax(sizezmax)<-s && max(newZmax)<M-s;       %Type 1: Exits bottom and doesn't hit top beforehand (increase L)
    Zmaxtype=1;
elseif newZmax(sizezmax)<=s && max(newZmax)>=M-s;  %Type 2: Hits or exits bottom, and hits top beforehand (keep L)
    Zmaxtype=2;
elseif newZmax(sizezmax)>s && min(newZmax)>s;      %Type 3: Exits top and doesn't hit bottom beforehand (decrease L)
    Zmaxtype=3;
elseif newZmax(sizezmax)>s && min(newZmax)<=s      %Type 4: Exits top and hits bottom beforehand (keep L)
    Zmaxtype=4;
elseif sizezmax==n+1-tau && tbottommax<sizezmax; %Type 4: Feasible but hits the bottom before the end time (keep L)
    Zmaxtype=4;
else
    Zmaxtype=5;                       %Either feasible or will need adjusting
end


% Defining Zmin
if newZmin(sizezmin)<-s && max(newZmin)<M-s;       %Type 1: Exits bottom and doesn't hit top beforehand (increase L)
    Zmintype=1;
elseif newZmin(sizezmin)<=s && max(newZmin)>=M-s;  %Type 2: Hits or exits bottom, and hits top beforehand (keep L)
    Zmintype=2;
elseif newZmin(sizezmin)>s && min(newZmin)>s;      %Type 3: Exits top and doesn't hit bottom beforehand (decrease L)
    Zmintype=3;
elseif newZmin(sizezmin)>s && min(newZmin)<=s      %Type 4: Exits top and hits bottom beforehand (keep L)
    Zmintype=4;
elseif sizezmin==n+1-tau && tbottommin<sizezmin; %Type 4: Feasible but hits the bottom before the end time (keep L)
    Zmintype=4;
else
    Zmintype=5;                       %Either feasible or will need adjusting
end


% Calculating Ztype (note we may not have covered all cases here, needs work)
if (Zmaxtype == 1) && (Zmintype == 1)   %Both exit bottom and doesn't hit top before (Increase L)
    Ztype = 1;
elseif (Zmaxtype == 3) && (Zmintype == 3)   %Both exit top and doesn't hit bottom before (decrease L)
    Ztype = 3;
elseif Zmaxtype==2 || Zmaxtype==4 || Zmaxtype==5    %Zmax is feasible (Keep L)
    Ztype = 4;
elseif Zmintype==2 || Zmintype==4 || Zmintype==5    %Zmin is feasible (Keep L)
    Ztype = 5;
else
    Ztype = 2; %Different constraints are broken (Keep L)
end

end


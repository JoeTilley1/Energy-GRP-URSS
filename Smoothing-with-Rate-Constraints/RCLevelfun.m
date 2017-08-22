function [Z_max,Z_min] = RCLevelfun(p,tau,X0,L,M,a1,a2,s,q1,q2)

%This program gives the level of the store according to the multiplier L from the end of time tau.  
%The last entry of Z is the first time that a capacity constraint is broken.
%Z does not include X0 (e.g. when tau=1, then X0 corresponds to time 0).
%% Defining variables

n=numel(p);
X=X0;   
Z=[];           %list storing level of store
x=[];           %list storing values of x
y=[];           %list storing values of y
Lvec(1) = L;    %lambda_1

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
    X=X+x(k);
    Z=[Z,X];
end

% % See the 'particle through a cave' analogy
% hold on
% plot(y)
% plot(alpha(q1,a1,a2).*q1-p)
% plot(alpha(-q2,a1,a2).*-q2-p)
% hold off


% We now modify y's unaffecting TV maximally and minimally
[y_min,y_max] = RCyminmax(y,p,tau,a1,a2,q1,q2);


% We calculate the associated x_max and x_min values to y_max and y_min
% (NOTE: To save computational time, this step should be incorporated into
% RCyminmax.m since we calculate these values at some point during that
% function. Could save up to 1/2 the running time.)

x_max = [];
x_min = [];

for i=1:length(y)
    if p(i+tau-1)+y_max(i)>=0
        x_max(i) = a1*(p(i+tau-1)+y_max(i));
    else
        x_max(i) = (p(i+tau-1)+y_max(i))/a2;
    end
    if p(i+tau-1)+y_min(i)>=0
        x_min(i) = a1*(p(i+tau-1)+y_min(i));
    else
        x_min(i) = (p(i+tau-1)+y_min(i))/a2;
    end
end

% % See picture
% hold on
% plot(y)
% plot(y_max)
% plot(y_min)
% plot((alpha(q1,a1,a2).*q1-p(tau:length(p))))
% plot((alpha(-q2,a1,a2).*-q2-p(tau:length(p))))

% Finding levels Z_max and Z_min according to these values of x_max and x_min

X_max=X0;   
Z_max=[];
for k=1:length(x_max)
    X_max=X_max+x_max(k);
    Z_max=[Z_max,X_max];
end

X_min=X0;   
Z_min=[];
for k=1:length(x_min)
    X_min=X_min+x_min(k);
    Z_min=[Z_min,X_min];
end

% % See the store level associated to y,ymax,ymin
% hold on
% plot(Z(1:max(length(Z_max),length(Z_min))))
% plot(Z_max)
% plot(Z_min)
end




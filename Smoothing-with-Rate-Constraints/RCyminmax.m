function [ymin,ymax] = RCyminmax(y,p,tau,a1,a2,q1,q2)
%RCYMINMAX: Given y, we find ymin and ymax such that the rate constraints
%are not broken and the TV is unchanged

%% Defining variables

%y has errors of ~e-16 where values should be equal, so we get rid of these
for i=2:length(y)
    if abs(y(i)-y(i-1))<10^(-10)
        y(i)=y(i-1);
    end
end

%Maximum and minimum attainable values of y
maxheight = alpha(q1,a1,a2)*q1-p(tau:length(p));
minheight = alpha(-q2,a1,a2)*-q2-p(tau:length(p));


% See the 'particle a cave' analogy
% % hold on
% % plot(y)
% % plot(alpha(q1,a1,a2).*q1-p(tau:length(p)))
% % plot(alpha(-q2,a1,a2).*-q2-p(tau:length(p)))


% Classifying minimum and maximum's of y
% For flat peaks/troughs, we store only the first value of the peak/trough

minimas = [];
maximas = [];

for i=2:length(y)-1
    initiali = i;
    if y(i)-y(i-1)>0 && y(i)-y(i+1)>=0
        while abs(y(i)-y(i+1))==0
            i=i+1;
            if i==length(y)
                maximas = [maximas, initiali];
                break
            end
        end
        if i~=length(y)
            if y(initiali)-y(i+1)>0
                maximas = [maximas, initiali];
            end
        end
    end
    if y(i)-y(i-1)<0 && y(i)-y(i+1)<=0
        while abs(y(i)-y(i+1))==0
            i=i+1;
            if i==length(y)
                minimas = [minimas, initiali];
                break
            end
        end
        if i~=length(y)
            if y(initiali)-y(i+1)<0
                minimas = [minimas, initiali];
            end
        end
    end
end

%Adding first value (if a list is empty, there is one extremum and it is at index=1)
if ~isempty(maximas) && ~isempty(minimas)
    if minimas(1)<maximas(1)
        maximas = [1,maximas];
    else
        minimas = [1,minimas];
    end
elseif isempty(maximas) && isempty(minimas)
    maximas = 1;    %Straight line case. WLOG we choose maximas = 1 instead of minimas
elseif isempty(maximas) && ~isempty(minimas)
    maximas = 1;
elseif ~isempty(maximas) && isempty(minimas)
    minimas = 1;
else
    disp('Error calculating ymin and ymax (classifying extremums)')
end

%Extremums
extr = sort([maximas, minimas]);

%Defining ymax,ymin
ymax = zeros(size(y));
ymin = ymax;
ymax(1) = y(1);
ymin(1) = y(1);

for i=2:length(extr)
    if ismember(extr(i),maximas)
        ymin(extr(i-1):extr(i)) = y(extr(i-1):extr(i));
        % Start of flip method
        % Reseting variables
        flippedp = 0; 
        flippedymax=0;
        Lvec= 0;
        x=0;
        flippedp = fliplr(p(extr(i-1)+tau-1:extr(i)+tau-1));
        flippedymax=0;
        flippedymax(1)= y(extr(i));
        Lvec(1) = 0; %This value is unused intentionally
        x(1) = 0; %This value is unused intentionally
        Lvec(2) = flippedymax(1);
        for k=2:length(flippedp)
            alphaxx = flippedp(k)+Lvec(k);
            if alphaxx <= -a2*q2
                x(k) = -q2;
            elseif alphaxx >= q1/a1
                x(k) = q1;
            elseif flippedp(k)+Lvec(k)>=0
                x(k) = a1*(flippedp(k)+Lvec(k));
            else
                x(k) = (flippedp(k)+Lvec(k))/a2;
            end
            flippedymax(k) = alpha(x(k),a1,a2)*x(k)-flippedp(k);
            Lvec(k+1) = flippedymax(k);
        end
        ymax(extr(i-1):extr(i)) = fliplr(flippedymax);
        %End of flip method   
    elseif ismember(extr(i),minimas)
        ymax(extr(i-1):extr(i)) = y(extr(i-1):extr(i));
        %Start of flip method
        %Reseting variables
        flippedp = 0;
        flippedymin=0;
        Lvec= 0;
        x=0;
        flippedp = fliplr(p(extr(i-1)+tau-1:extr(i)+tau-1));
        flippedymin(1)= y(extr(i));
        Lvec(1) = 0; %This value is unused intentionally
        x(1) = 0; %This value is unused intentionally
        Lvec(2) = flippedymin(1);
        for k=2:length(flippedp)
            alphaxx = flippedp(k)+Lvec(k);
            if alphaxx <= -a2*q2
                x(k) = -q2;
            elseif alphaxx >= q1/a1
                x(k) = q1;
            elseif flippedp(k)+Lvec(k)>=0
                x(k) = a1*(flippedp(k)+Lvec(k));
            else
                x(k) = (flippedp(k)+Lvec(k))/a2;
            end
            flippedymin(k) = alpha(x(k),a1,a2)*x(k)-flippedp(k);
            Lvec(k+1) = flippedymin(k);
        end
        ymin(extr(i-1):extr(i)) = fliplr(flippedymin);
        %End of flip method
    else
        disp('Error calculating ymin/ymax')
    end  
end

%Making first values all equal
ymax(1) = y(1);
ymin(1) = y(1);

%After last peak/trough, y=ymin=ymax
ymin(extr(length(extr)):length(y))=y(extr(length(extr)):length(y));
ymax(extr(length(extr)):length(y))=y(extr(length(extr)):length(y));

% % See picture
% hold on
% plot(y)
% plot(ymax)
% plot(ymin)
% plot((alpha(q1,a1,a2).*q1-p(tau:length(p))))
% plot((alpha(-q2,a1,a2).*-q2-p(tau:length(p))))

end


function output = alpha(x,a1,a2)
%ALPHA: This function calculates the coefficient of efficiency
    if x>=0;
        output=1/a1;
    else
        output=a2;
    end
end


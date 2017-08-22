function Z=JoeLevelfun(p,tau,X0,L,M,a1,a2,s)

%This program gives the level of the store according to the multiplier L from the end of time tau.  
%The last entry of Z is the first time that a capacity constraint is broken.
%Z does not include X0 (e.g. when tau=1, then X0 corresponds to time 0).

n=numel(p);
X=X0;   
Z=[];

%Since the properties of the gradient of a stretched string match those of the usual multipler conditions, we have L=alpha(x_t)*x_t-w_t.
for i=1:n+1-tau
    pii=p(i+tau-1);     
    if L+pii>=0;
        x=(L+pii)*a1;
    else
        x=(L+pii)/a2;
    end
    X=X+x;
    Z=[Z,X];
    if X<-s || X>M+s;
        break
    end
end

end


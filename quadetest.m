function quadetest(x,varargin)
%QUADETEST: Quade test for non parametric two way ANalysis Of VAriance.
%This function performs the Quade test to analyze unreplicated complete block
%designs.
%Dana Quade in 1979 proposed a test that is often more powerful than the
%Friedman test. It also eliminates block differences but weights the raw data
%indicate possibly more marked treatment effects. Whereas the Friedman test is
%basically an extension of the sign test, the Quade test is effectively an
%extension of the Wilcoxon signed rank test and is equivalent to it when the
%treatments are two.
%
% Syntax: 	STATS=quadetest(X,alpha)
%      
%     Inputs:
%           X - data matrix
%           ALPHA - significance level (default = 0.05).
%     Outputs:
%           Quade Statistic
%           Multiple comparisons (eventually)
%
%      Example: 
%
% x=[115 142 36 91 28; 28 31 7 21 6; 220 311 108 51 117; 82 56 24 46 33;...
% 256 298 124 46 84; 294 322 176 54 86; 98 87 55 84 25];
%
%           Calling on Matlab the function: quadetest(x)
%
%           Answer is:
%
% QUADE TEST FOR IDENTICAL TREATMENT EFFECTS:
% TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS
% --------------------------------------------------------------------------------
% Number of observation: 35
% Number of blocks: 7
% Number of treatments: 5
% --------------------------------------------------------------------------------
% F-statistic approximation
% Quade test statistic W: 10.3788
% F=W df-num=4 df-denom=24 - p-value (2 tailed): 0.0001
% --------------------------------------------------------------------------------
%  
% POST-HOC MULTIPLE COMPARISONS
% --------------------------------------------------------------------------------
% Critical value: 35.6981
% Absolute difference among mean ranks
%      0     0     0     0     0
%     18     0     0     0     0
%     51    69     0     0     0
%     69    87    18     0     0
%     63    81    12     6     0
% 
% Absolute difference > Critical Value
%      0     0     0     0     0
%      0     0     0     0     0
%      1     1     0     0     0
%      1     1     0     0     0
%      1     1     0     0     0
%
%           Created by Giuseppe Cardillo
%           giuseppe.    cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). QUADETEST: Quade test for non parametric two way ANalysis Of VAriance
% http://www.mathworks.com/matlabcentral/fileexchange/25926

%Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty'}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p,x,varargin{:});
assert(all(x(:,2) == fix(x(:,2))),'Warning: all elements of column 2 of input matrix must be whole numbers')
alpha=p.Results.alpha;
clear p

[r,c]=size(x); %dimension of the input matrix
R=zeros(r,c); %preallocation
%For each block, compute the ranks
for I=1:r
    R(I,:)=tiedrank(x(I,:));
end
%Compute the range of each block and then rank them.
Q=tiedrank(range(x,2));
%Compute a modified version of the Friedman matrix
rij=(R-(c+1)/2).*repmat(Q,1,c);
Ti=sum(rij);
T2=sum(Ti.^2);
rij2=sum(sum(rij.^2));
T3=T2/r;
T4=rij2-T3;
k=r-1;
W=k*T3/T4; %The Quade statistic.
%The Quade statistic is approximable with the F distribution.
dfn=c-1;
dfd=dfn*k;
p=1-fcdf(W,dfn,dfd);

%display results
tr=repmat('-',1,80); %set the divisor
disp('QUADE TEST FOR IDENTICAL TREATMENT EFFECTS: TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS')
disp(tr)
disp(table(r*c,r,c,'VariableNames',{'Observations','Blocks','Treatments'}))
disp('QUADE''S STATISTICS: F-statistic approximation')
disp(tr)
disp(table(W,dfn,dfd,p,'VariableNames',{'W','DF_numerator','DF_denominator','two_tailed_p_value'}))
if p<alpha
    disp(' ')
    disp('POST-HOC MULTIPLE COMPARISONS')
    disp(tr)
    tmp=repmat(Ti,c,1); Rdiff=abs(tmp-tmp'); %Generate a matrix with the absolute differences among ranks
    cv=tinv(1-alpha/2,dfd)*realsqrt(2*r*T4/dfd); %critical value
    mc=Rdiff>cv; %Find differences greater than critical value
    %display results
    fprintf('Critical value: %0.4f\n',cv)
    disp('Absolute difference among mean ranks')
    %disp(tril(Rdiff))
    disp(tril(Rdiff))
    disp('Absolute difference > Critical Value')
    disp(tril(mc))
end

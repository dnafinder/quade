function [stats, mc] = quadetest(x, varargin)
%QUADETEST Quade test for nonparametric two-way ANOVA (complete block design).
%
%   [STATS, MC] = QUADETEST(X)
%   [STATS, MC] = QUADETEST(X, 'Alpha', ALPHA, 'PostHoc', POSTHOC, 'Display', DISPLAY)
%
%   This function performs the Quade test to analyze unreplicated complete
%   block designs.
%
%   Dana Quade (1979) proposed a test that is often more powerful than the
%   Friedman test. It also eliminates block differences but weights the
%   raw data to emphasise blocks indicating more marked treatment effects.
%   Whereas the Friedman test is basically an extension of the sign test,
%   the Quade test is effectively an extension of the Wilcoxon signed rank
%   test and is equivalent to it when there are only two treatments.
%
%   Inputs:
%     X        - data matrix (Blocks x Treatments), real, finite, non-NaN.
%                Each row is a block, each column a treatment.
%
%   Name–Value pair arguments:
%     'Alpha'  - significance level (scalar in (0,1), default = 0.05).
%
%     'PostHoc'- logical-like flag to enable post-hoc multiple comparisons:
%                true  -> perform multiple comparisons when global H0 is rejected
%                false -> only perform the global Quade test
%                (default = true)
%
%     'Display'- logical-like flag controlling command-window output:
%                true  -> print tables and messages (default)
%                false -> run silently, only return outputs
%
%   Outputs:
%     STATS - structure with test results:
%       .nObs          - number of observations (blocks * treatments)
%       .blocks        - number of blocks (rows)
%       .treatments    - number of treatments (columns)
%       .W             - Quade statistic
%       .F_df_num      - numerator degrees of freedom
%       .F_df_denom    - denominator degrees of freedom
%       .F_p           - two-tailed p-value for F approximation
%       .alpha         - significance level
%       .rejectNull    - true if the treatments do not have identical effects
%
%     MC   - structure with post-hoc multiple comparison results (empty if
%            PostHoc is false or the global null is not rejected):
%       .method      - 'Quade-Conover-type LSD'
%       .Rdiff       - matrix of absolute differences among treatment scores
%       .cv          - critical value for differences
%       .pvalue      - [] (no exact p-values computed)
%       .significant - logical matrix (lower triangle) indicating significant
%                      pairwise differences at the chosen alpha
%
%   Example:
%
%     x = [115 142  36  91  28;
%           28  31   7  21   6;
%          220 311 108  51 117;
%           82  56  24  46  33;
%          256 298 124  46  84;
%          294 322 176  54  86;
%           98  87  55  84  25];
%
%     [stats, mc] = quadetest(x);
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, this would be an appropriate format:
%   Cardillo G. (2009). QUADETEST: Quade test for non parametric two way
%   ANalysis Of VAriance. Available on GitHub:
%   https://github.com/dnafinder/quade

% -------------------------------------------------------------------------
% Input Error handling
% -------------------------------------------------------------------------
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'x', @(z) validateattributes(z, {'numeric'}, ...
    {'2d','real','finite','nonnan','nonempty'}, mfilename, 'X', 1));
addParameter(p, 'Alpha',   0.05, @(z) validateattributes(z, {'numeric'}, ...
    {'scalar','real','finite','nonnan','>',0,'<',1}, mfilename, 'Alpha'));
addParameter(p, 'PostHoc', true,  @validateLogicalLike);
addParameter(p, 'Display', true,  @validateLogicalLike);
parse(p, x, varargin{:});

x          = p.Results.x;
alpha      = p.Results.Alpha;
postHocFlg = logical(normalizeLogicalLike(p.Results.PostHoc));
displayFlg = logical(normalizeLogicalLike(p.Results.Display));

% -------------------------------------------------------------------------
% Basic dimensions
% -------------------------------------------------------------------------
[r, c] = size(x);            % r = blocks, c = treatments
nObs   = r * c;

% -------------------------------------------------------------------------
% Quade transformation
% -------------------------------------------------------------------------
R = zeros(r, c);             % ranks within each block
for ii = 1:r
    R(ii, :) = tiedrank(x(ii, :));
end

% Rank the block ranges
Q = tiedrank(range(x, 2));   % weight for each block (column vector, r x 1)

% Modified Friedman-type matrix
rij  = (R - (c + 1) / 2) .* repmat(Q, 1, c);  % weighted centred ranks
Ti   = sum(rij);                              % treatment scores (1 x c)
T2   = sum(Ti.^2);
rij2 = sum(sum(rij.^2));
T3   = T2 / r;
T4   = rij2 - T3;
k    = r - 1;

% Quade statistic (F-approximation)
W   = k * T3 / T4;
dfn = c - 1;
dfd = dfn * k;
pval = 1 - fcdf(W, dfn, dfd);

% -------------------------------------------------------------------------
% Build STATS structure
% -------------------------------------------------------------------------
stats = struct();
stats.nObs       = nObs;
stats.blocks     = r;
stats.treatments = c;
stats.W          = W;
stats.F_df_num   = dfn;
stats.F_df_denom = dfd;
stats.F_p        = pval;
stats.alpha      = alpha;
stats.rejectNull = (pval < alpha);

% -------------------------------------------------------------------------
% Display results (if requested)
% -------------------------------------------------------------------------
mc = struct([]);

if displayFlg
    tr = repmat('-', 1, 80);
    disp('QUADE TEST FOR IDENTICAL TREATMENT EFFECTS: TWO-WAY BALANCED, COMPLETE BLOCK DESIGNS');
    disp(tr);
    disp(table(nObs, r, c, 'VariableNames', {'Observations','Blocks','Treatments'}));
    disp('QUADE''S STATISTICS: F-statistic approximation');
    disp(tr);
    disp(table(W, dfn, dfd, pval, ...
        'VariableNames', {'W','DF_numerator','DF_denominator','two_tailed_p_value'}));
    
    if ~stats.rejectNull
        % Nessun commento aggiuntivo nell’originale, ma potresti aggiungerlo se vuoi
    end
end

% -------------------------------------------------------------------------
% Post-hoc multiple comparisons (if requested and global H0 rejected)
% -------------------------------------------------------------------------
if stats.rejectNull && postHocFlg
    if displayFlg
        disp(' ');
        disp('POST-HOC MULTIPLE COMPARISONS');
        tr = repmat('-', 1, 80);
        disp(tr);
    end
    
    % Generate matrix of absolute differences among treatment scores
    tmp   = repmat(Ti, c, 1);
    Rdiff = abs(tmp - tmp');
    
    % Critical value (Quade-Conover-type LSD)
    cv     = tinv(1 - alpha/2, dfd) * realsqrt(2 * r * T4 / dfd);
    mcMask = Rdiff > cv;
    
    if displayFlg
        fprintf('Critical value: %0.4f\n', cv);
        disp('Absolute difference among mean ranks');
        disp(tril(Rdiff));
        disp('Absolute difference > Critical Value');
        disp(tril(mcMask));
    end
    
    mc(1).method      = 'Quade-Conover-type LSD';
    mc(1).Rdiff       = Rdiff;
    mc(1).cv          = cv;
    mc(1).pvalue      = [];
    mc(1).significant = tril(mcMask, -1);
end

% -------------------------------------------------------------------------
% Output behaviour
% -------------------------------------------------------------------------
if nargout == 0
    clear stats mc
end

end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------
function tf = validateLogicalLike(x)
%VALIDATELOGICALLIKE Helper for inputParser: check logical-like values.
    try
        normalizeLogicalLike(x);
        tf = true;
    catch
        tf = false;
    end
end

function y = normalizeLogicalLike(x)
%NORMALIZELOGICALLIKE Convert various logical-like inputs to true/false.
    if islogical(x)
        y = x;
    elseif isnumeric(x) && isscalar(x)
        y = (x ~= 0);
    elseif ischar(x) || (isstring(x) && isscalar(x))
        s = lower(char(x));
        if any(strcmp(s, {'true','on','yes'}))
            y = true;
        elseif any(strcmp(s, {'false','off','no'}))
            y = false;
        else
            error('Invalid logical-like value: %s', s);
        end
    else
        error('Invalid type for logical-like option.');
    end
end

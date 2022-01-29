function normalized_weight = robust_weight(x)
% This function returns normalized feature weights based on robust
% regression. It takes one input argument that is the expression data (as numeric matrix), the
% rows contain feature and the columns contain samples. The values in the
% expression data should not be negative.

mean_f = mean(x,2);
cv_f = std(x,[],2)./mean_f;
log2_mean = log2(mean_f);
log2_cv = log2(cv_f);
robust_fit = robustfit(log2_mean,log2_cv);

%% assigning feature weights
robust_fit = round(robust_fit*100)/100;
cv_fit = 2.^(robust_fit(1)+(log2_mean * robust_fit(2)));
log2_cv_fit = log2(cv_fit);
weight = log2_cv - log2_cv_fit;


%% weights normalization
contains_neg = any(weight<0);
if contains_neg == 1
     min_abs_weight = abs(min(weight));
     pos_weight = weight + min_abs_weight;
     normalized_weight = pos_weight/(sum(pos_weight));
else
    normalized_weight = weight/(sum(weight));
end
end

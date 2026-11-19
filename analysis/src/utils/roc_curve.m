function [false_positive_rate, true_positive_rate, thresholds, auc] = roc_curve(labels, scores)
%ROC_CURVE Calculate a binary ROC curve for scores where larger is positive.

labels = logical(labels(:));
scores = scores(:);
valid = isfinite(scores);
labels = labels(valid);
scores = scores(valid);

n_positive = sum(labels);
n_negative = sum(~labels);
if n_positive == 0 || n_negative == 0
    false_positive_rate = [0; 1];
    true_positive_rate = [0; 1];
    thresholds = [Inf; -Inf];
    auc = NaN;
    return
end

[sorted_scores, order] = sort(scores, 'descend');
sorted_labels = labels(order);
last_in_tie = [find(diff(sorted_scores) ~= 0); numel(sorted_scores)];
true_positive = cumsum(sorted_labels);
false_positive = cumsum(~sorted_labels);

true_positive_rate = [0; true_positive(last_in_tie) ./ n_positive];
false_positive_rate = [0; false_positive(last_in_tie) ./ n_negative];
thresholds = [Inf; sorted_scores(last_in_tie)];
auc = trapz(false_positive_rate, true_positive_rate);
end

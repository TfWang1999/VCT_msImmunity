function metrics = fit_metrics(observed, predicted)
%FIT_METRICS Calculate common goodness-of-fit summaries.

observed = observed(:);
predicted = predicted(:);
valid = isfinite(observed) & isfinite(predicted);
observed = observed(valid);
predicted = predicted(valid);

residual = observed - predicted;
metrics.n = numel(observed);
metrics.rmse = sqrt(mean(residual.^2));
metrics.mae = mean(abs(residual));
metrics.nrmse = metrics.rmse / max(range(observed), eps);
ss_res = sum(residual.^2);
ss_tot = sum((observed - mean(observed)).^2);
if ss_tot > 0
    metrics.r2 = 1 - ss_res / ss_tot;
else
    metrics.r2 = NaN;
end
end


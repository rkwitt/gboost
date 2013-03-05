function [classifier, cfun] = gboost1d5 (G, Y, nu, conv_epsilon);

% 1.5-class graph based LPBoosting

disp(['=== GLPBOOST = BEGIN = 1.5-class graph based LPBoosting ===']);
% FIXME: make convergence threshold (0.1) configurable
[classifier, cfun] = lpboost (G, Y, conv_epsilon, nu, ...
	@(X, Y, u, beta, max_col) findhypothesis_graph (X, Y, u, beta, max_col, 1), 1);
disp(['=== END ===']);


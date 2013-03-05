function [classifier, cfun] = gboost2 (G, Y, nu, conv_epsilon);

% 2-class graph based LPBoosting

disp(['=== GLPBOOST = BEGIN = 2-class graph based LPBoosting ===']);
% FIXME: make convergence threshold (0.1) configurable
[classifier, cfun] = lpboost (G, Y, conv_epsilon, nu, ...
	@(X, Y, u, beta, max_col) findhypothesis_graph (X, Y, u, beta, max_col, 2), 2);
disp(['=== END ===']);



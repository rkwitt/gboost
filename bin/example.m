% Example training and testing a 2-class graph boosting classifier

disp(['Loading example graphs...']);
load example-graphs.mat

disp(['   ', num2str(length(train_G)), ' training samples']);
disp(['   ', num2str(length(test_G)), ' test samples']);
disp(' ');

disp(['We train a 2-class graph boosting classifier.']);
disp(['the settings are:']);
disp(['   nu = 0.2, the LPBoost nu-parameter controlling training accuracy']);
disp(['   conv_epsilon = 0.05, the LPBoost convergence tolerance parameter']);
disp(['   max_col = 25, generate 25 hypotheses in each iteration (multiple pricing)']);
disp(' ');
disp(['Please press return to start training...']);
pause;

[cl, cfun] = gboost2 (train_G, train_Y, 0.2, 0.05);
disp(['The classifier has been trained successfully.']);
disp(['There are ', num2str(length(find(cl.alpha > 1e-5))), ' active subgraph stumps.']);
disp(['Please press return to test the classifier...']);
pause;

[Yout, Yreal, GY] = cfun (test_G);
accuracy = length(find(Yout == test_Y))/length(test_Y);
disp(['   test accuracy = ', num2str(accuracy)]);
[auc, eer, curve] = rocscore (Yreal, test_Y);
disp(['   test ROC AUC = ', num2str(auc)]);
disp(['   test ROC EER = ', num2str(eer)]);


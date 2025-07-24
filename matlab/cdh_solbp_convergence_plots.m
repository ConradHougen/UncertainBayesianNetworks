% File for post-processing convergence of variances vs means
clear;
linf = 5000;

diffpso = zeros(linf, 6);
diffvso = zeros(linf, 6);

load('.\Post_Fusion\sdiam_nodes1^7_iter50.mat');
diffpso(:,1) = (infpso - infpsolbp);
diffvso(:,1) = (infvso - infvsolbp);

load('.\Post_Fusion\sdiam_nodes1^7_iter100.mat');
diffpso(:,2) = (infpso - infpsolbp);
diffvso(:,2) = (infvso - infvsolbp);

load('.\Post_Fusion\sdiam_nodes1^7_iter150.mat');
diffpso(:,3) = (infpso - infpsolbp);
diffvso(:,3) = (infvso - infvsolbp);

load('.\Post_Fusion\sdiam_nodes1^7_iter200.mat');
diffpso(:,4) = (infpso - infpsolbp);
diffvso(:,4) = (infvso - infvsolbp);

load('.\Post_Fusion\sdiam_nodes1^7_iter250.mat');
diffpso(:,5)= (infpso - infpsolbp);
diffvso(:,5) = (infvso - infvsolbp);

load('.\Post_Fusion\sdiam_nodes1^7_iter300.mat');
diffpso(:,6) = (infpso - infpsolbp);
diffvso(:,6) = (infvso - infvsolbp);

% load('.\Post_Fusion\sdiam_nodes1^7_iter35.mat');
% diffpso(:,7) = (infpso - infpsolbp);
% diffvso(:,7) = (infvso - infvsolbp);
% 
% load('.\Post_Fusion\sdiam_nodes1^7_iter40.mat');
% diffpso(:,8) = (infpso - infpsolbp);
% diffvso(:,8) = (infvso - infvsolbp);
% 
% load('.\Post_Fusion\sdiam_nodes1^7_iter45.mat');
% diffpso(:,9) = (infpso - infpsolbp);
% diffvso(:,9) = (infvso - infvsolbp);
% 
% load('.\Post_Fusion\sdiam_nodes1^7_iter50.mat');
% diffpso(:,10) = (infpso - infpsolbp);
% diffvso(:,10) = (infvso - infvsolbp);

dpso_pp = diffpso;
%dpso_pp = rmoutliers(diffpso);
figure;
boxplot(dpso_pp, 'Labels', {'50','100','150','200','250','300'});
title('Error vs. Num Iterations (Means)');

dvso_pp = diffvso;
%dvso_pp = rmoutliers(diffvso);
figure;
boxplot(dvso_pp, 'Labels', {'50','100','150','200','250','300'});
title('Error vs. Num Iterations (Vars)');

fcorrect_pso = sum(diffpso < 1e-8, 1)/linf;
fcorrect_vso = sum(diffvso < 1e-8, 1)/linf;
figure; plot(fcorrect_pso, 'o'); title('Fraction Good Mean Inferences vs. Num Iterations');
figure; plot(fcorrect_vso, 'o'); title('Fraction Good Var Inferences vs. Num Iterations');




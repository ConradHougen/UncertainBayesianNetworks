% post processing script to check if "outlier" inferences cause decbod
% error
clear;
load('.\Post_Fusion\21nodenet_iter30.mat');

figure;
scatter(infpso, infpsolbp, 'MarkerEdgeColor','b');
title('Inferred Means: 21 node iter30');
xlabel('Mean Value (SPN Method)');
ylabel('Mean Value (SOLBP)');

figure;
scatter(infvso, infvsolbp, 'MarkerEdgeColor','b');
title('Inferred Variances: 21 node iter30');
xlabel('Variance Value (SPN Method)');
ylabel('Variance Value (SOLBP)');


infpsodiff = infpso-infpsolbp;
figure;
histogram(infpsodiff);
title('Inferred Means Error Distribution');
xlabel('Error');
ylabel('Frequency');

infvsodiff = infvso-infvsolbp;
figure;
histogram(infvsodiff);
title('Inferred Variances Error Distribution');
xlabel('Error');
ylabel('Frequency');


% Now remove center lobe of mean inferences and plot decbod
outer = find(abs(infpsodiff) > 0.001);
nouter = length(outer);
pct_outer = nouter/length(infpso) * 100;
fprintf(1, "Using %1.2f pct of total inferences for decbod plot\n", pct_outer);

weso_outer = cdh_gen_subjective_opinions(infpso(outer), infvso(outer));
wesolbp_outer = cdh_gen_subjective_opinions(infpsolbp(outer), infvsolbp(outer));
[pa3, pp3] = perfbound_beta_new2(pgt(outer), weso_outer);
[pa4, pp4] = perfbound_beta_new2(pgt(outer), wesolbp_outer);

figure; 
plot(pp3, pa3, 'b', pp4, pa4, 'r+', pp4, pp4, 'k--');
title('DeCBoD Outer Only: 21 node iter30');
xlabel('Desired Confidence');
ylabel('Actual Confidence');
legend('SPN Method', 'SOLBP', 'Reference Line', 'Location', 'northwest');

figure; 
plot(pp, pa, 'b', pp2, pa2, 'r+', pp2, pp2, 'k--');
title('DeCBoD Normal: 21 node iter30');
xlabel('Desired Confidence');
ylabel('Actual Confidence');
legend('SPN Method', 'SOLBP', 'Reference Line', 'Location', 'northwest');
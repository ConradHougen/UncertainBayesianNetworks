% Script to run sequence of experiments and save results

%=======================================================
clear;
FRACOBS = 0.10;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_10pctrandom.mat');

clear;
FRACOBS = 0.10;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_10pctrandom.mat');

clear;
FRACOBS = 0.10;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_10pctrandom.mat');

clear;
FRACOBS = 0.10;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_10pctrandom.mat');

clear;
FRACOBS = 0.10;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_10pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.2;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_20pctrandom.mat');

clear;
FRACOBS = 0.2;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_20pctrandom.mat');

clear;
FRACOBS = 0.2;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_20pctrandom.mat');

clear;
FRACOBS = 0.2;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_20pctrandom.mat');

clear;
FRACOBS = 0.2;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_20pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.30;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_30pctrandom.mat');

clear;
FRACOBS = 0.30;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_30pctrandom.mat');

clear;
FRACOBS = 0.30;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_30pctrandom.mat');

clear;
FRACOBS = 0.30;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_30pctrandom.mat');

clear;
FRACOBS = 0.30;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_30pctrandom.mat');
%=======================================================\
%=======================================================
clear;
FRACOBS = 0.40;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_40pctrandom.mat');

clear;
FRACOBS = 0.40;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_40pctrandom.mat');

clear;
FRACOBS = 0.40;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_40pctrandom.mat');

clear;
FRACOBS = 0.40;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_40pctrandom.mat');

clear;
FRACOBS = 0.40;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_40pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.50;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_50pctrandom.mat');

clear;
FRACOBS = 0.50;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_50pctrandom.mat');

clear;
FRACOBS = 0.50;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_50pctrandom.mat');

clear;
FRACOBS = 0.50;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_50pctrandom.mat');

clear;
FRACOBS = 0.50;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_50pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.60;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_60pctrandom.mat');

clear;
FRACOBS = 0.60;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_60pctrandom.mat');

clear;
FRACOBS = 0.60;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_60pctrandom.mat');

clear;
FRACOBS = 0.60;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_60pctrandom.mat');

clear;
FRACOBS = 0.60;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_60pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.70;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_70pctrandom.mat');

clear;
FRACOBS = 0.70;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_70pctrandom.mat');

clear;
FRACOBS = 0.70;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_70pctrandom.mat');

clear;
FRACOBS = 0.70;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_70pctrandom.mat');

clear;
FRACOBS = 0.70;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_70pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.80;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_80pctrandom.mat');

clear;
FRACOBS = 0.80;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_80pctrandom.mat');

clear;
FRACOBS = 0.80;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_80pctrandom.mat');

clear;
FRACOBS = 0.80;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_80pctrandom.mat');

clear;
FRACOBS = 0.80;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_80pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 0.90;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_90pctrandom.mat');

clear;
FRACOBS = 0.90;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_90pctrandom.mat');

clear;
FRACOBS = 0.90;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_90pctrandom.mat');

clear;
FRACOBS = 0.90;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_90pctrandom.mat');

clear;
FRACOBS = 0.90;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_90pctrandom.mat');
%=======================================================
%=======================================================
clear;
FRACOBS = 1;
netfile = '3node_1_1_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_chain_100pctrandom.mat');

clear;
FRACOBS = 1;
netfile = '3node_1_2.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tree_100pctrandom.mat');

clear;
FRACOBS = 1;
netfile = '3node_2_1.file';
cdh_run_decbod_experiments;
save('Experiments/3n_V_100pctrandom.mat');

clear;
FRACOBS = 1;
netfile = 'net2_dag.file';
cdh_run_decbod_experiments;
save('Experiments/9n_dag_100pctrandom.mat');

clear;
FRACOBS = 1;
netfile = '3node_triangle.file';
cdh_run_decbod_experiments;
save('Experiments/3n_tri_100pctrandom.mat');
%=======================================================
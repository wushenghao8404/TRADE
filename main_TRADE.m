%% Calling the solvers
% For large population sizes, consider using the Parallel Computing Toolbox of MATLAB.
% Else, the program can be slow.
function [  ] = main_TRADE(benchmark_name,prob_index,reps)
    addpath('.\BasicFunction\');
	if strcmp(benchmark_name,'GECCO20MaTO')
		Tasks = GECCO20MaTO_benchmark(prob_index);
		ntask = length(Tasks);
		% parameter setting for the benchmark
		ps= 100 * ntask; % population size for multitasking  
		change_stage_gen = 100;
		ncheckpts = 1000;
		maxfes = ntask*1e5;
	elseif  strcmp(benchmark_name,'CEC19MaTO')
		Tasks = CEC19MaTO_benchmark(prob_index);
		ntask = length(Tasks);
		% parameter setting for the benchmark
		ps= 100 * ntask; % population size for multitasking  
		change_stage_gen = 100;
		ncheckpts = 1000;
		maxfes = ntask * 1e5; 
	end
	
	disp(['benchmark ' benchmark_name ' problem ' num2str(prob_index)]);
	data_MTDE = TRADE( Tasks,ps,maxfes,ncheckpts,change_stage_gen,reps );
end


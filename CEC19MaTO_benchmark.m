function [Tasks, g] = benchmark(index)
%BENCHMARK function
%   Input
%   - index: the index number of problem set
%
%   Output:
%   - Tasks: benchmark problem set

    nTask=50;
    switch(index)                      
        case 1 
            load('.\CEC19ManyTasks\GoTask1.mat');  % loading data from folder .\Tasks
            load('.\CEC19ManyTasks\RotationTask1.mat');
            dim = 50;
            
            for i=1:nTask
                Tasks(i).dim = dim;   
                Tasks(i).fnc = @(x)Rosenbrock(x,RotationTask1{i},GoTask1(i,:));
                Tasks(i).Lb=-50*ones(1,dim);   
                Tasks(i).Ub=50*ones(1,dim);   
                g(i,:)=-1*ones(1,dim)+GoTask1(i,:);
            end                        
        case 2
            load('.\CEC19ManyTasks\GoTask2.mat');  % loading data from folder .\Tasks
            load('.\CEC19ManyTasks\RotationTask2.mat');
            dim = 50;
            
            for i=1:nTask
                Tasks(i).dim = dim;   
                Tasks(i).fnc = @(x)Ackley(x,RotationTask2{i},GoTask2(i,:));
                Tasks(i).Lb=-50*ones(1,dim);   
                Tasks(i).Ub=50*ones(1,dim);   
                g(i,:)=GoTask2(i,:);
            end   
        case 3
            load('.\CEC19ManyTasks\GoTask3.mat');  % loading data from folder .\Tasks
            load('.\CEC19ManyTasks\RotationTask3.mat');
            dim = 50;
            
            for i=1:nTask
                Tasks(i).dim = dim;   
                Tasks(i).fnc = @(x)Rastrigin(x,RotationTask3{i},GoTask3(i,:));
                Tasks(i).Lb=-50*ones(1,dim);   
                Tasks(i).Ub=50*ones(1,dim);   
                g(i,:)=GoTask3(i,:);
            end               
        case 4
            load('.\CEC19ManyTasks\GoTask4.mat');  % loading data from folder .\Tasks
            load('.\CEC19ManyTasks\RotationTask4.mat');
            dim = 50;
            
            for i=1:nTask
                Tasks(i).dim = dim;   
                Tasks(i).fnc = @(x)Griewank(x,RotationTask4{i},GoTask4(i,:));
                Tasks(i).Lb=-100*ones(1,dim);   
                Tasks(i).Ub=100*ones(1,dim);   
                g(i,:)=GoTask4(i,:);
            end   
        case 5
            load('.\CEC19ManyTasks\GoTask5.mat');  % loading data from folder .\Tasks
            load('.\CEC19ManyTasks\RotationTask5.mat');
            dim = 50;
            
            for i=1:nTask
                Tasks(i).dim = dim;   
                Tasks(i).fnc = @(x)Weierstrass(x,RotationTask5{i},GoTask5(i,:));
                Tasks(i).Lb=-0.5*ones(1,dim);   
                Tasks(i).Ub=0.5*ones(1,dim);   
                g(i,:)=GoTask5(i,:);
            end   
        case 6
            load('.\CEC19ManyTasks\GoTask6.mat');  % loading data from folder .\Tasks
            load('.\CEC19ManyTasks\RotationTask6.mat');
            dim = 50;
            
            for i=1:nTask
                Tasks(i).dim = dim;   
                Tasks(i).fnc = @(x)Schwefel(x,RotationTask6{i},GoTask6(i,:));
                Tasks(i).Lb=-500*ones(1,dim);   
                Tasks(i).Ub=500*ones(1,dim);   
                g(i,:)=-420.9687*ones(1,dim)+GoTask6(i,:);
            end               
            
    end
end
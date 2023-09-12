function [ ng,gId ] = task_grouping( evbf_tasks )
%UNTITLED2 此处显示有关此函数的摘要
    [ntask,~]=size(evbf_tasks);
    task_param=[];
    for taskid=1:ntask
        evbf=evbf_tasks(taskid,:);
        lb=mean(evbf)-4*std(evbf);
        evbf=max(1e-25,evbf-lb);
        evbf=log(evbf);
        task_param=[task_param ;evbf];
    end
    [ ng ] = ddCRP( task_param );
    [gId] = kmeans(task_param,ng);
	gId=gId';
end


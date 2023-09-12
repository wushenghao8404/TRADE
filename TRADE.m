function [ data_MTDE ] = TRADE( Tasks,ps,maxfes,ncheckpts,change_stage_gen,nrun )
udim=max([Tasks.dim]);
ntask=length(Tasks);
lb=0;
ub=1;
EvBestFitness=1e25*ones(nrun*ntask,ncheckpts);
GroupId=zeros(nrun,ntask);
Gbx=ones(nrun*ntask,udim);
EvParamId=zeros(nrun,ntask);
EvParamCount={};
EvSucParamCount={};
RunTime=zeros(nrun,1);

for run=1:nrun
	tic
	disp(['run ' num2str(run)]);
	poolF = [0.5];
	poolFSize = length(poolF);
	poolCr = [0.1 0.5 0.9];
	poolCrSize = length(poolCr);
    with_ar = true;
    subps = ceil(ps/ntask)*ones(1,ntask);
    lbF = 0.01;
    ubF = 1.0;
    lbCr = 0.01;
    ubCr = 1.0;
    mF = 0.5*ones(1,ntask); % sensitive setting
    mCr = 0.5*ones(1,ntask); % sensitive setting (e.g., on Schwefel, (0.95,0.05) works better than (0.5,0.5))
    sdF = 0.1*ones(1,ntask);
    sdCr = 0.1*ones(1,ntask);
    gbfit = 1e25*ones(1,ntask);
    gbx = ones(ntask,udim);
	
	% count the used times of different base solvers
	paramId = zeros(1,ntask);  % selected base solver ID
	paramCount = zeros(3,ntask);
	sucParamCount = zeros(3,ntask);
	EvParamCount{run}=zeros(ntask * 3,ncheckpts);
	EvSucParamCount{run}=zeros(ntask * 3,ncheckpts);
    
    tp = 0.05;
    gen = 1;
    fes = 0;
    evQuality = zeros(1,ntask);  % evolution quality
	stage = 1;
    pseet = 0.0;
    
	% initialization of populations
    for i=1:ntask
        lpop(i).pop = rand(subps(i),udim)*(ub-lb)+lb;
        lpop(i).fit = zeros(1,subps(i));
        for j=1:subps(i)
            lpop(i).fit(j) = fnceval(Tasks(i),lpop(i).pop(j,:));
            fes=fes+1;
        end
        
        arc(i).pop = lpop(i).pop;
        arc(i).fit = lpop(i).fit;
        [gbfit(i),gbi] = min(lpop(i).fit);
        gbx(i,:) = lpop(i).pop(gbi,:);
        
        EvBestFitness(ntask*(run-1)+i,1)=gbfit(i);
        
        sortpop(i).pop=[];
    end
    
    while fes<maxfes

	   if stage == 2
			pseet = 0.5 + 0.5 * (gen-change_stage_gen) / (1000-change_stage_gen);
	   end
       disp(['TRADE gen ' num2str(gen) ' fes ' num2str(fes) ' avgfit ' num2str(mean(gbfit)) ' stage ' num2str(stage)]);
	   gen=gen+1;
       
       for i=1:ntask
           % evolution of single task
           narx=size(arc(i).pop,1);
           if narx>subps(i)
                index=randperm(narx);
                arc(i).pop=arc(i).pop(index(1:subps(i)),:);
                arc(i).fit=arc(i).fit(index(1:subps(i)));
           end
           
           % sort the population from best to worst
           [~,index]=sort(lpop(i).fit);
           sortpop(i).pop=lpop(i).pop(index,:);
           
           if with_ar
                pop_arx(i).pop=[lpop(i).pop; arc(i).pop];
            else
                pop_arx(i).pop=lpop(i).pop;
            end
       end
       
       for i=1:ntask
           vpop = zeros(subps(i),udim);
           upop = zeros(subps(i),udim);
           ufit = zeros(1,subps(i));
		   
		   if stage == 2
			    gid = groupId(i);
			    gs = groupsize(gid);
                if gs>poolCrSize
	                good_taskids = groupMemberId{gid}(task_rank(groupMemberId{gid})<gs/poolCrSize);
                else
  	                good_taskids = groupMemberId{gid};
                end
				% good_taskids = groupMember{gid}(task_rank(groupMember{gid})<task_rank(i));
			    good_tasks_size = length(good_taskids);
		   end
		   
           for j=1:subps(i)
				if stage==1
					% stage 1 - use JADE as the base solver with the associated parameter setting to collect information data for task grouping
					F=randn()*sdF(i)+mF(i);
                    F=min(ubF,max(lbF,F));
					Cr=randn()*sdCr(i)+mCr(i);
                    Cr=min(ubCr,max(lbCr,Cr));
				else
					% stage 2 - combining self evolution mechanism and transferring successful parameters from other tasks
					F=randn()*sdF(i)+mF(i);
					F=min(ubF,max(lbF,F));
					
					% generate Cr for i-th individual
					if rand()<pseet && rand()>1/task_rank(i) && groupsize(gid)>1
						% learn from successful evolution experience of other tasks
						seltaskid = good_taskids(randi(good_tasks_size));
						[~,selParamId] = max(sucParamCount(:,seltaskid)./(paramCount(:,seltaskid)+1e-10));
						paramCount(selParamId, i) = paramCount(selParamId, i) + 1;
						
	                    Cr=randn()*0.1+poolCr(selParamId);
						Cr=min(ubCr,max(lbCr,Cr));
					else
						selParamId = paramId(i);
						paramCount(selParamId, i) = paramCount(selParamId, i) + 1;
						
						Cr=randn()*sdCr(i)+mCr(i);
						Cr=min(ubCr,max(lbCr,Cr));
					end
				end

                r1=randi(subps(i));
                while r1==j
                    r1=randi(subps(i));
                end
                
                npop_arx=size(pop_arx(i).pop,1);
                r2=randi(npop_arx);
                while r2==j||r2==r1
                    r2=randi(npop_arx);
                end

                % mutation
				pbi=randi(ceil(tp*subps(i)));
				vpop(j,:)=lpop(i).pop(j,:)+F*(sortpop(i).pop(pbi,:)-lpop(i).pop(j,:))...
						 +F*(lpop(i).pop(r1,:)-pop_arx(i).pop(r2,:));

                % crossover
                drnd = randi(udim);
                for d=1:udim
                    if d==drnd|| rand()<Cr
                        upop(j,d)=vpop(j,d);
                    else
                        upop(j,d)=lpop(i).pop(j,d);
                    end
                end

                % constrain handling
                upop(j,:)=min(max(upop(j,:),lb),ub);

                % evaluation
                ufit(j)=fnceval(Tasks(i),upop(j,:));
                fes=fes+1;

                % selection
                if ufit(j)<lpop(i).fit(j)
                    arc(i).pop=[arc(i).pop; lpop(i).pop(j,:)];
                    arc(i).fit=[arc(i).fit lpop(i).fit(j)];
                    lpop(i).pop(j,:)=upop(j,:);
                    lpop(i).fit(j)=ufit(j);
					
					if stage == 2
						sucParamCount(selParamId, i) = sucParamCount(selParamId, i) + 1;
					end
                    if ufit(j)<gbfit(i)
                        gbfit(i)=ufit(j);
                        gbx(i,:)=upop(j,:);
                    end
                end 
           end
       end
       
       % record mF, mCr, gbest, number of used different parameter settings
       for ii=1:ntask
           EvBestFitness(ntask*(run-1)+ii,gen)=gbfit(ii);
       end
	   EvSucParamCount{run}(:,gen)=reshape(sucParamCount,3*ntask,1);
	   EvParamCount{run}(:,gen)=reshape(paramCount,3*ntask,1);
       
	   % evolution quality evaluation
	   if stage == 2
		   for ii=1:ntask
			   diff1 = abs(EvBestFitness(ntask*(run-1)+ii,change_stage_gen)-EvBestFitness(ntask*(run-1)+ii,gen));
			   diff2 = abs(EvBestFitness(ntask*(run-1)+ii,1)-EvBestFitness(ntask*(run-1)+ii,gen));
			   evQuality(ii)=diff1/(diff2+1e-25);
		   end
	   end
       
        % switch stage
        if gen==change_stage_gen
            stage=2;
            task_rank=ones(1,ntask);
            [ ngroup,groupId ] = task_grouping( EvBestFitness(ntask*(run-1)+1:ntask*run,1:change_stage_gen) );
            groupsize = zeros(1,ngroup);
            groupMemberId = cell(1,ngroup);
            GroupId(run,:)=groupId;
            for gi=1:ngroup
                groupMemberId{gi} = find(groupId==gi);
                groupsize(gi)=length(groupMemberId{gi});
            end
            
            disp(['stage ' num2str(stage) ' number of groups ' num2str(ngroup)]);
            for gi=1:ngroup
                disp(groupMemberId{gi});
            end
            
            for i=1:ntask
                mF(i)=poolF(randi(poolFSize));
                paramId(i)=randi(poolCrSize);
                mCr(i)=poolCr(paramId(i));
            end
            EvParamId(run,:)=paramId;
        end
	   
	    % update task rank
        if stage==2
            for gi=1:ngroup
               [~,ind]=sort(-evQuality(groupMemberId{gi}));
               [~,task_rank(groupMemberId{gi})]=sort(ind);
            end
        end
    end
    used_time=toc;
    RunTime(run)=used_time;
    disp(['Running time: ' num2str(used_time)]);
    disp('Gbest fitness: ');
    disp(gbfit);
    Gbx((run-1)*ntask+1:(run-1)*ntask+ntask,:)=gbx;
end
data_MTDE.EvBestFitness=EvBestFitness;
data_MTDE.GroupId=GroupId;
data_MTDE.Gbx=Gbx;
data_MTDE.EvSucParamCount=EvSucParamCount;
data_MTDE.EvParamCount=EvParamCount;
data_MTDE.EvParamId=EvParamId;
data_MTDE.RunTime=RunTime;
end


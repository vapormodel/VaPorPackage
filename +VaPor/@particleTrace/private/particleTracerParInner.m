function [T_TransientStore] = particleTracerParInner(...
    particlesToInsert,InletPoints,InletFlows, OutletPoints, NumGM_WM,...
    NumArt, Tree, TotalTime, Timestep, SaveTime, InjectionLength, Option_MinHops, minHops)
%PARTICLETRACERPARINNER Summary of this function goes here
%   Detailed explanation goes here
% Create the vectors to hold the inputs. 
    times_taken = zeros(ceil(particlesToInsert/10000),1);
    particles_generated = times_taken;
    limit_counter = 0;
start_cycle = tic;
T_TransientStore = zeros(size(Tree,1),TotalTime/SaveTime);
%fid = fopen('hello.txt', 'a+');
for npnts = 1:particlesToInsert
    if mod(npnts, 100000) == 0
%         time_to_now = toc(start_cycle);
%         time_per_comped_point = time_to_now/npnts;
%         remaining_points = particlesToInsert - npnts;
%         est_time_to_comp = remaining_points*time_per_comped_point/60;
limit_counter = limit_counter+1;
particles_generated(limit_counter) = npnts;
times_taken(limit_counter) = toc(start_cycle);

fit_res = polyfit(particles_generated(1:limit_counter),times_taken(1:limit_counter),2);

total_exp_time = polyval(fit_res, particlesToInsert);

est_time_to_comp = total_exp_time-times_taken(limit_counter);
        
        % Report results
        fprintf('Completed %d of %d points (%.2f %%) [TE: %.2f min ETR: %.2f min]\n', npnts, particlesToInsert, (npnts/particlesToInsert)*100, times_taken(limit_counter)/60, est_time_to_comp/60);
    elseif mod(npnts, 10000) == 0
        % Calculate the Estimated Time Remaining:
        %         time_to_now = toc(start_cycle);
        %         time_per_comped_point = time_to_now/npnts;
        %         remaining_points = particlesToInsert - npnts;
        %         est_time_to_comp = remaining_points*time_per_comped_point/60;
        
        % Report results
        limit_counter = limit_counter+1;
particles_generated(limit_counter) = npnts;
times_taken(limit_counter) = toc(start_cycle);
        fprintf('Completed %d of %d points (%.2f %%)\n', npnts, particlesToInsert, (npnts/particlesToInsert)*100);
    end
    
    AllocatedRows = 1000;
    travelLog = zeros(AllocatedRows,2);
    rowcount=1;
    %tic
    % point path function
    % choose inlet
    if length(InletPoints) ~= 1
        location = InletPoints(find(cumsum(InletFlows/sum(InletFlows))>=rand,1,'first'))+NumGM_WM;
    else
        location = InletPoints + NumGM_WM;
    end
    %location = InletPoints(1)+NumGM_WM;
    
    travelLog(rowcount,:) = [location(1),0];
    rowcount=rowcount+1;
    % update location
    
    % Minhops Parameters:
    enteredPorous = false;
    numberOfHops = 0;
    
    while  ~ismember(location,OutletPoints+NumArt+NumGM_WM)
        location_matrix = Tree{location};
        if isempty(location_matrix)
            break
        end
        location_new = location_matrix(find(cumsum(location_matrix(:,2))>=rand,1,'first'),[1,3]);
        
        % If the particle hasn't passed through the a number of voxels then
        % it should be prevented from leaving.
        if Option_MinHops
            % Check if we have entered the porous media for the first time.
            if ~enteredPorous && (location_new(1,1) <= NumGM_WM)
                % We have entered the porous media for the first time.
                enteredPorous = true;
                numberOfHops = numberOfHops + 1;
            elseif enteredPorous && (numberOfHops < minHops)
                % We havent taken the minimum number of steps.
                loopCheck = 0;
                escapeLoop = false;
                while ~escapeLoop
                    location_new = location_matrix(find(cumsum(location_matrix(:,2))>=rand,1,'first'),[1,3]);
                    loopCheck = loopCheck + 1;
                    if loopCheck > 10000
                        %error('Upper Ceiling on Iterations Reached');
                        escapeLoop = true;
                        enteredPorous = false;
                    end
                    if location_new(1,1) <= NumGM_WM
                        escapeLoop = true;
                    end
                    if size(location_matrix, 1) == 1
                        % exit the loop
                        escapeLoop = true;
                        enteredPorous = false;
                    end
                    % Check that there is the possibility of transiting
                end
                numberOfHops = numberOfHops + 1;
            end
            
        end
        
        %             sum_col4 = sum(location_matrix(:,4));
        %             if sum_col4==0, sum_col4=1; end
        %             [max_val,max_ind] = max(location_matrix(:,2)-location_matrix(:,4)/sum_col4);
        %             location_new = location_matrix(max_ind,[1,3]);
        %             Tree{location}(max_ind,4) = Tree{location}(max_ind,4) + 1;
        if size(location_new,1) == 0
            break
        end
        travelLog(rowcount,:) = [location_new(1,1),0];
        
        % update time
        if location_new(1,1) <= NumGM_WM
            %                         fudged_MTT = location_new(1,2) + ((2*rand)-1)*location_new(1,2);
            fudged_MTT = normrnd(location_new(1,2),abs(location_new(1,2)/4));
        else
            fudged_MTT = location_new(1,2);
        end
        %         fudged_MTT = location_new(1,2);
        if  location_new(2) < 0
            travelLog(rowcount-1,2) = travelLog(rowcount-1,2) + abs(fudged_MTT);
        else
            travelLog(rowcount,2) = abs(fudged_MTT);
        end
        
        %             travelLog(rowcount,2) = travelLog(rowcount,2) + ((2*rand)-1)*travelLog(rowcount,2);
        
        location = location_new(1,1);
        rowcount=rowcount+1;
    end
    
    if rowcount > AllocatedRows
        disp(['TravelLog larger than allocated rows ( ' num2str(rowcount) ' > ' num2str(AllocatedRows) ')'])
    end
    
    travelLog = travelLog(logical(travelLog(:,1)),:);
    
    
    %tlstring = sprintf('%d,' , travelLog(:,1));
    %tlstring = tlstring(1:end-1);% strip final comma
    %fprintf(fid, '%s\n', tlstring);
    
    
    
    timelist(npnts,1) = sum(travelLog(:,2));
    timelist(npnts,2) = sum(travelLog(1:find(travelLog(:,1)<NumGM_WM,1,'first')-1,2));
    timelist(npnts,3) = sum(travelLog(find(travelLog(:,1)<NumGM_WM,1,'first'):find(travelLog(:,1)<=NumGM_WM,1,'last'),2));
    timelist(npnts,4) = sum(travelLog(find(travelLog(:,1)<=(NumGM_WM),1,'last')+1:end,2));
    
    cumm_time = cumsum(travelLog(:,2));
    %         con_domains(travelLog(1,1),1) = con_domains(travelLog(1,1),1) + 1;
    
    for N = 1:(TotalTime/Timestep)
        t = N*Timestep;
        
        if t-InjectionLength <= cumm_time(end)
            
            %             bar1 = find(cumm_time>=t,1,'first');
            %             con_domains(travelLog(bar1,1),n+1) = con_domains(travelLog(bar1,1),n+1) + 1;
            
            bar1 = find(cumm_time>=t,1,'first'); if isempty(bar1), bar1 = numel(cumm_time); end
            bar2 = find(cumm_time>=(t-InjectionLength),1,'first');
            %             con_domains(travelLog(bar1:bar2,1),n+1) = con_domains(travelLog(bar1:bar2,1),n+1) + 1;
            
            if bar1 == bar2
                if bar1 == numel(cumm_time)
                    T_TransientStore(travelLog(bar1,1),N) = T_TransientStore(travelLog(bar1,1),N) + min((cumm_time(end)-(t-InjectionLength)),InjectionLength);
                else
                    T_TransientStore(travelLog(bar1,1),N) = T_TransientStore(travelLog(bar1,1),N) + min(t,InjectionLength);
                end
            else
                T_TransientStore(travelLog(bar1,1),N) = T_TransientStore(travelLog(bar1,1),N) + t-cumm_time(bar1-1);
                T_TransientStore(travelLog(bar2,1),N) = T_TransientStore(travelLog(bar2,1),N) + cumm_time(bar2)-max((t-InjectionLength),0);
                T_TransientStore(travelLog(bar2+1:bar1-1,1),N) = T_TransientStore(travelLog(bar2+1:bar1-1,1),N) + cumm_time(bar2+1:bar1-1)-cumm_time(bar2:bar1-2);
            end
            
        end
    end
end
time_to_now = toc(start_cycle)/60;
% Report results
fprintf('Completed %d points in %.2f min\n', particlesToInsert, time_to_now);
%fclose(fid);
end


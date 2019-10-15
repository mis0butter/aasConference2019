% Monte carlos 

% Add current folder to path 
dir_code = pwd; 
addpath(dir_code); 

% # of test cases 
n_cases = 1000;  

% Create overall summary 
date_folder = sprintf('outputs/%s_MC_summary.txt', datestr(now, 'yyyy-mmm-dd_HH.MM.SS')); 
MC_summary = fopen(date_folder, 'w'); 

% run tests 
for i = 1:n_cases 
    
    fprintf('Case %d \n', i); 
    
    % Main test script. Contains: 
    % - inputs.m
    % - steering_profile.m
    % - post_processing.m
    main
    
    % Create folder to put results 
    cd('outputs'); 
    date_str = datestr(now, 'yyyy-mmm-dd_HH.MM.SS'); 
    mkdir(date_str); 
    cd(date_str); 
    
    % Create summary.txt 
    create_summary
    
    % Save figures and plots 
    savefig('profile_unit_sphere'); 
    plots = allchild(0); 
    for k = 1:length(plots)
        saveas(plots(k), strcat(get(plots(k), 'Name'), '.png'))
    end 
    
    % Save workspace 
    close all;
    save workspace.mat 
    
    % Return to code folder 
	cd(dir_code) 
    
    % Write into MC_summary
    err_final = acosd(dot(Pf_N, P_phi3_N(end, :)));
    fprintf(MC_summary, '%s \t \t', date_str); 
    fprintf(MC_summary, 'Final error: %.2f deg \n', err_final); 

end 

fclose(MC_summary); 
disp('All done!') 

%% Plot summary of data 

date_folder = '2019-Oct-06_12.41.38_MC_summary.txt'; 
date_folder = sprintf('outputs/%s', date_folder); 

% Import data
MC_data = importdata(date_folder); 

% Save only error values 
for i = 1:length(MC_data)
    cell_contents = MC_data{i}; 
    MC_data{i} = cell_contents(38:41); 
end 
MC_data = cellfun(@str2double, MC_data); 

% Plot 
figure
    scatter( 1:1:length(MC_data), MC_data, 'filled')  
    grid on; 
    ylabel('deg') 
    xlabel('Test Case') 
    title('Final Error') 

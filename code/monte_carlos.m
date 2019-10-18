% Monte carlos 

% Add current folder to path 
dir_code = pwd; 
addpath(dir_code); 

% # of test cases 
n_cases = 3;  

% Create overall summary 
folder_date_string = datestr(now, 'yyyy-mmm-dd_HH.MM.SS'); 
date_folder = sprintf('outputs/%s_MC_summary.txt', folder_date_string); 
MC_summary = fopen(date_folder, 'w'); 

% Create large errors summary 
errors_folder = sprintf('outputs/%s_large_errors.txt', folder_date_string); 
errors_summary = fopen(errors_folder, 'w'); 
fprintf(errors_summary, '\t \t \t \t \t \t '); 
fprintf(errors_summary, 'Pi_N \t \t \t \t \t \t Pf_N \t \t \t \t \t \t S_N \t \t \t \t \t \t '); 
fprintf(errors_summary, 'alpha \t \t \t a_max \t \t \t w_max \t \t \t '); 
fprintf(errors_summary, 'phi_1 \t \t \t phi_2 \t \t \t error \n'); 


% run tests 
for i = 1:n_cases 
    
    fprintf('Case %d \n', i); 
    
    % Main test script. Contains: 
    % - inputs.m
    % - steering_profile.m
    % - post_processing.m
    main()
    
    % Create folder to put results 
    cd('outputs'); 
    date_str = datestr(now, 'yyyy-mmm-dd_HH.MM.SS'); 
    mkdir(date_str); 
    cd(date_str); 
    
    % Create summary.txt 
    create_summary()
    
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
    
    % Write into errors_summary
    if err_final > 0
        fprintf(errors_summary, '%s \t ', date_str); 
        fprintf(errors_summary, '[%.3f, %.3f, %.3f] \t ', Pi_N(1), Pi_N(2), Pi_N(3)); 
        fprintf(errors_summary, '[%.3f, %.3f, %.3f] \t ', Pf_N(1), Pf_N(2), Pf_N(3)); 
        fprintf(errors_summary, '[%.3f, %.3f, %.3f] \t ', S_N(1), S_N(2), S_N(3)); 
        fprintf(errors_summary, '%.3f rad \t ', alpha); 
        fprintf(errors_summary, '%.3f rad/s^2 \t ', aMax); 
        fprintf(errors_summary, '%.3f rad/s \t ', wMax); 
        fprintf(errors_summary, '%.3f rad \t \t ', phi1); 
        fprintf(errors_summary, '%.3f rad \t \t ', phi2); 
        fprintf(errors_summary, '%.3f rad \t \t ', phi3); 
        fprintf(errors_summary, '%.3f deg \n', err_final); 
    end 

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

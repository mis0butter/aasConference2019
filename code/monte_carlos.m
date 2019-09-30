% Monte carlos 

% Add current folder to path 
dir_code = pwd; 
addpath(dir_code); 

% # of test cases 
n = 100; 

% Create overall summary 
MC_summary = fopen(sprintf('outputs/%s_MC_summary.txt', datestr(now, 'yyyy-mmm-dd_HH.MM.SS')), 'w'); 

% run tests 
for i = 1:n 
    
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


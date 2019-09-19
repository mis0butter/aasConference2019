% Monte carlos 

% Add current folder to path 
addpath('.'); 

% # of test cases 
n = 5; 

% Create overall summary 
MC_summary = fopen(sprintf('outputs/%s_MC_summary.txt', datestr(now, 'dd-mmm-yyyy_HH.MM.SS')), 'w'); 

% run tests 
for i = 1:n 
    
    % Main test script 
    main
    
    % Create folder to put results 
    cd outputs; 
    date_str = datestr(now, 'dd-mmm-yyyy_HH.MM.SS'); 
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
	cd ../..
    
    % Write into MC_summary
    err_final = acosd(dot(Pf_N, P_phi3_N(end, :)));
    fprintf(MC_summary, '%s \t \t', date_str); 
    fprintf(MC_summary, 'Final error: %.2f deg \n \n', err_final); 

    
end 


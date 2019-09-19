% Monte carlos 

% Add current folder to path 
addpath('.'); 

% # of test cases 
n = 5; 

% run tests 
for i = 1:n 
    
    % Main test script 
    main
    
    % Create folder to put results 
    cd outputs; 
    dir_output = datestr(now, 'dd-mmm-yyyy_HH.MM.SS'); 
    mkdir(dir_output); 
    cd(dir_output); 
    
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
    
	cd ../..
    
end 


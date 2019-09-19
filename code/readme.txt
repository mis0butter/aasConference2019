SAS Algorithm Code 

Run main.m for stand-alone simulation 

Run run_monte_carlos.m to run batch tests 
    also runs create_summary 
    saves plots, workspace vars, and profile unit sphere fig to outputs folder 

To toggle plot generation: 
    in post_processing.m: 
        set plot_option = 1 or 0 for various plots 

Command window output: 
    phi 1+2+3 
        sum of slew angles 
    Final error 
        Diff of target Pf and propagated final attitude 

To set alpha = 0 
    in sun_vector.m: 
        un-comment "Make alpha = 0" section 

Workspace for results in paper: 
    alpha0.mat 
    alphaNot0.mat 
SAS Algorithm Code 

Run main.m 

To generate plots: 
    in post_processing.m: 
        set plot_option = 1 for various plots 

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
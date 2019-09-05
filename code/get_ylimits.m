function ylimits = get_ylimits(data) 
% Set y limits of axis 

% rng = range(data(:)); 

rng = max(data(:)) - min(data(:)); 

midp = min(data(:)) + rng/2; 
ylimits = [ midp - 1.2*rng/2, midp + 1.2*rng/2 ]; 

end 
function brownian_motion = brownian_motion(dim, times)
    % dim: how many dimensions?
    % times: sorted row vector (1xT) of positive timestamps to sample the bm
    % :returns borwnian_motion: (dim x T) vector where every column i is the
    %   bm at time times(i)
    if isa(times, 'cell')
        [flat_times, ~, indices] = unique(cat(2, times{:}));
        lengths = cellfun(@length, times);
    else
        flat_times = times;
    end
    
    s = size(flat_times);
    std_normal = randn(dim, s(2)); % s(2) = T (1:columsize, 2:rowsize)
    sqrt_diffs = sqrt(diff(cat(2,0, flat_times))); % insert t_0 = 0
    flat_brownian_motion = cumsum(std_normal .* sqrt_diffs, 2);
    
    if isa(times, 'cell')
        non_unique = flat_brownian_motion(:,indices);
        slices = cat(1, 0, cumsum(lengths));
        brownian_motion = cell(length(lengths),1);
        for i = 1:(length(slices)-1)
           brownian_motion{i} = non_unique(:,(slices(i)+1):slices(i+1));
        end
    else
        brownian_motion = flat_brownian_motion;
    end
 end

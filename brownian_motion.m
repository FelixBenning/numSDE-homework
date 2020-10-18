function brownian_motion = brownian_motion(dim, times)
    % dim: how many dimensions?
    % times: sorted row vector (1xT) of positive timestamps to sample the bm
    % :returns borwnian_motion: (dim x T) vector where every column i is the
    %   bm at time times(i)
    s = size(times);
    std_normal = randn(dim, s(2)); % s(2) = T (1:columsize, 2:rowsize)
    sqrt_diffs = sqrt(diff(cat(2,0, times))); % insert t_0 = 0
    brownian_motion = cumsum(std_normal .* sqrt_diffs, 2);
 end

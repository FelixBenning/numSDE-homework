function conf_interval = confidence_interval(avg, variance, niveau, N)
    conf_interval = avg+norminv([niveau/2, 1-niveau/2])*sqrt(variance/N);
end
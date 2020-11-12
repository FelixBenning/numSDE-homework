%% variable definitions

fun = @(x) x-x^3;
x_0 = 0;
v_0 = 0;
eta = 0.6;
sigma = 0.2;
T = 100;

%% (i)

deltas = arrayfun( @(n) 10^-n, -1:2);

elapsed = zeros(1,length(deltas));
x = cell(1, length(deltas));
v = cell(1, length(deltas));

for idx = 1:length(deltas)
    tic
    [x{idx}, v{idx}] = leapfrog(fun, deltas(idx), x_0, v_0, T, eta, sigma);
    elapsed(idx) = toc;
end

%% delta = 10
plot(x{1}, v{1});

%% delta = 1
plot(x{2}, v{2});

%% delta = 0.01
plot(x{4}, v{4});

%% delta = 0.1
plot(x{3}, v{3});

%% (ii)

N = 10^3;
delta = 1;

v_avg = zeros(1, N);
x_neg1_1 = zeros(1,N);
for idx = 1:N
    [x,v] = leapfrog(fun, delta, x_0, v_0, T, eta, sigma);
    v_avg(idx) = mean(v);
    x_big_enough = x(x>-1);
    x_satisfied = x_big_enough(x_big_enough<1);
    x_neg1_1(idx) = length(x_satisfied)/length(x);
end

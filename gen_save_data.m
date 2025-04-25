% Generating data object for saving the results.

data.pdf = fs_idle;
data.cdf = Fs_idle;
data.buffer_levels = xs;
data.config = config;

matrices.Qs = Qs;
matrices.Q1s = Q1s;
matrices.Q10 = Q10;
matrices.Rs = Rs;
matrices.R1s = R1s;
data.matrices = matrices;

KPIs.delay_avg = mean_buffer/config.params.rate;
KPIs.power_avg = avg_power;
KPIs.time_act = time_act;
KPIs.time_deact = time_deact;

data.KPIs = KPIs;
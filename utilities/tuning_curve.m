function out = tuning_curve(spk, pos, bins, dt)

pos_spike = pos(spk>0);
sp_count = histcounts(pos_spike,bins);
time_count = dt*histcounts(pos,bins);

w = gausswin(15);
w = w/sum(w);
fr = sp_count./rep_zero(time_count);
fr_sm = conv(fr,w,'same');

out = fr_sm;
end
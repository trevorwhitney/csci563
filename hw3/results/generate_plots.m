close all;
colors = ['r', 'g', 'b'];
scripts = {'broadcast'; 'broadcast_mpi'; 'broadcast_pipe'};
for i=0:2
    str = sprintf('%s_results.csv', char(scripts(i+1)));
    s = dlmread(str, ',', 1, 0);
    x = s(:,1);
    y = s(:,2);
    plot(x,y,colors(i+1));
    hold on;
end

xlabel('# of Processors');
ylabel('Time (seconds)');
title('Comparing Implementations of MPI Broadcast');

legend('Binomial Brodcast', 'MPI\_Bcast()', 'Pipeline Broadcast', 'Location', 'East');
close all;
colors = ['r', 'g', 'b'];
for i=0:2
    str = sprintf('%d_results.csv', i);
    s = dlmread(str, ',', 1, 0);
    x = s(:,1);
    y = s(:,2);
    plot(x,y,colors(i+1));
    hold on;
end

xlabel('# of Processors');
ylabel('Time (seconds)');
title('Sieve of Eratosthenes Results');

legend('Original Script', 'No Evens', 'No Evens, No Prime Comm.', 'Location', 'NorthWest');
function plot_different_time_steps(data1, data2, data3, time1, time2, time_step)
    
    a = time1/time_step;
    b = time2/time_step;

    n = abs(b - a);
    t = linspace(time1, time2, n);
    plot(t, data1((a + 1):b))
    hold on
    plot(t, data2((a + 1):b))
    hold on
    plot(t, data3((a + 1):b))
    legend("h = 0.01", "h= 0.005", "h = 0.001")
    ylabel("Probabilities.") % switch to voltage when comparing U
    xlabel("Time in milliseconds.")

end


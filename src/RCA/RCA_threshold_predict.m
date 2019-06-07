function snr_out = RCA_threshold_predict(B,snr_start,p,R,snr_R)

snr_step_array = [10, 1, 0.1, 0.01];
snr = snr_start;

for istep = 1:length(snr_step_array)
    snr_step = snr_step_array(istep);
    go = 1;
    while(go)
        snr = snr - snr_step;
        %����
        [flag] = RCA_predict(B,snr,p,R,snr_R,22);
        if flag == 0 
            go = 0;
            snr_out = snr + snr_step;
        end
    end
    snr = snr + snr_step;
end
end


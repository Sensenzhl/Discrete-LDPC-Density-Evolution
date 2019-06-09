matrix = [];
delta = 0.01;
for a = -38:0.1:38
    for b = -38:0.1:38
        R_original = 2 * atanh(tanh(a/2)*tanh(b/2));
        R = quantize_message(R_original,0.01);
        matrix = [matrix; [a,b,R]];
    end    
end


function [R]=quantize_message(R_original,delta)
    if R_original >= delta / 2
        R = delta * floor(R_original/delta + 1/2);
    elseif R_original <= -delta / 2
        R = delta * floor(R_original/delta - 1/2);
    else
        R = 0;
    end
end
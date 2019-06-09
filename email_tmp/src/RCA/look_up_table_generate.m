range = 40;
delta = 0.1;

p1 = ones(1, range/0.1 + 1);
p2 = ones(1, range/0.1 + 1);
a = zeros(1, range/0.1 + 1);

k_reg = [];
kmin = 0;
kmax = 0;

i_idx = 1;
j_idx = 1;

for i = 0 : delta : range
    for j = 0 : delta : range 
        i_tmp = i - range / 2;
        j_tmp = j - range / 2;      

        tmp = 2 * atanh(tanh(i_tmp/2 * delta) * tanh(j_tmp/2 * delta));
        if tmp >= delta / 2
            q = floor(tmp/delta + 0.5) * delta;
        elseif tmp <= -delta / 2
            q = ceil(tmp/delta - 0.5) * delta;
        else
            q = 0;
        end
               
        k = q / delta;
        
        if k < kmin
            kmin = k;
        end
            
        if k > kmax
            kmax = k;
        end
        
        k = k + 10;
        
        %k_idx = 1;
        
        a(k) = a(k) + p1(i_idx)*p2(j_idx);
        
        j_idx = j_idx + 1;
    end
    j_idx = 1;
    i_idx = i_idx + 1;
end
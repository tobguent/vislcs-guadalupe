function [C,A] = colormapping(F, MinValue, MaxValue, Name)

    if strcmp(Name, 'blue2red')
        t = min(max(0.0, (double(F) - MinValue) / (MaxValue - MinValue)), 1.0);
        R = (((5.0048.*t -8.0915).*t + 1.1657).*t + 1.4380).*t + 0.6639;
        G = (((7.4158.*t -15.9415).*t + 7.4696).*t + 1.2767).*t -0.0013;
        B = (((6.1246.*t -16.2287).*t + 11.9910).*t -1.4886).*t + 0.1685;
        A = abs(t - 0.5) * 2;
    end
    
    if strcmp(Name, 'white2blue')
        t = min(max(0.0, (double(F) - MinValue) / (MaxValue - MinValue)), 1.0);
        R = (((0.6599.*t + 0.5815).*t -1.7724).*t -0.4004).*t + 0.9636;
        G = (((0.7582.*t -1.3170).*t + 0.2448).*t -0.4792).*t + 0.9830;
        B = (((-1.6662.*t + 2.5921).*t -1.4702).*t -0.0314).*t + 0.9977;
        A = t;
    end
    
    if strcmp(Name, 'white2red')
        t = min(max(0.0, (double(F) - MinValue) / (MaxValue - MinValue)), 1.0);
        R = (((0.6786.*t -2.6266).*t + 1.6476).*t -0.2985).*t + 1.0031;
        G = (((-0.7348.*t + 3.0499).*t -2.9471).*t -0.3204).*t + 0.9594;    
        B = (((-3.0188.*t + 6.8674).*t -4.1791).*t -0.5610).*t + 0.9422;
        A = t;
    end
    
    if strcmp(Name, 'yellow')
        t = min(max(0.0, (double(F) - MinValue) / (MaxValue - MinValue)), 1.0);
        R = ones(size(t));
        G = ones(size(t));
        B = ones(size(t)) * 0.5;
        A = t;
    end
    
    C = cat(3,cat(3,R,G),B);
end


function s = detectGround(vt,a,epsi1,epsi2,T,k_curr)

% Calculation of noisy averages

vt_avg1 = mean(vt(k_curr-T-1:k_curr-1));
vt_avg2 = mean(vt(k_curr-T-a-1:k_curr-a-1));
vt_avg3 = mean(vt(k_curr-T-2*a-1:k_curr-2*a-1));

if abs(vt_avg3-vt_avg2)/a < epsi1 && abs(vt_avg1-vt_avg2)/a > epsi2

    s = 1; % ground change is detected
else 
    s = 0; % ground change is NOT detected


end
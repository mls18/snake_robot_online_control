function  [cn, ct ] = updateNominal(x_a,u_bar,vt,ct_tilde,cn_tilde,m,l,Ts,N)

% Parameters defining average vt
kao = (13.9*pi/180)^2*120*pi/180;
kdelta = 0;

% snake robot parameters
A = zeros(N-1, N);
D = zeros(N-1, N);

for i = 1: N-1
    
    A(i,i)= 1;
    A(i,i+1) = 1;
    
    D(i,i) = 1;
    D(i, i+1) = -1;
end
D_bar = D'*inv(D*D');
AD_bar = A*D_bar;
for i = 1 : N-1
    for j = 1 : N-1
        kdelta = kdelta + AD_bar(i,j)*sin((j-i)*40*pi/180);
    end
end
% Calculation of the average of vt and rho where ct = rho*cn

vt_avg_c = abs(mean(vt));
rho = (kao*kdelta)/(4*l*N*vt_avg_c+kao*kdelta);

% computation of cn using update formula
dv_phi = (x_a(2,N:end)-x_a(1,N:end))/Ts;

term1 = -x_a(1,N:end)'/m + ((1-rho)*vt(end-1)*A*D'*x_a(1,1:N-1)')/(2*l*m);
term2 = dv_phi'- u_bar - (cn_tilde/m*x_a(1,N:end)'+(ct_tilde-cn_tilde)/(2*l*m)*vt(end-1)*A*D'*x_a(1,1:N-1)');

cn = norm(term2)/norm(term1);
ct = rho*cn;

end

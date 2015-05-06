pkg load odepkg
gamma = 4;
omega = 2;
alpha = 2;
beta = 0;
% omega0 = 4.0;
dirname = 'data';

T = linspace(-5,5,501);
for jj=2:length(T)
	fprintf(1,'\rii = %d', jj);
	Opincs(:,:,jj-1) = propper2(T(jj), T(jj-1), gamma, omega, alpha, beta, omega0);

end
fprintf(1, '\r');
for ii=1:length(T);
	fprintf(1,'\rii = %d', ii);
	Op = eye(2,2);
	for jj=ii+1:length(T)
		Opinc = Opincs(:,:,jj-1);
		Op = Opinc*Op;
		Amp(ii,jj) = max(abs(eig(Op)));
	end
end
fprintf(1, '\r');
imagesc(T, T, log10(Amp'))
axis xy
colorbar
hold on;
[C, ~] = contour(T, T, log10(Amp'),[0, 0], 'k-', 'LineWidth', 2);
hold off;
drawnow;

fname = sprintf('Amp_omega0_%3.3f.mat', omega0);
save(fname, 'T', 'Amp', 'omega0');

cfname = sprintf('%s/C_omega0_%3.3f.mat', dirname, omega0);
cf = fopen(cfname, 'w');
while (~isempty(C))
        Nc = C(2,1);
        tmp = [C(:,2:Nc+1); omega0*ones(1,Nc)]; 
        fprintf(cf, '%f\t%f\t%f\n', tmp(:,2:end)); 
        fprintf(cf, '\n');
        C = C(:,Nc+2:end);
end
fclose(cf);


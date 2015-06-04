% Export magnitudes

fname = sprintf('data_%i_voltagem.data', figno);

myfile=fopen(fname,'w');
fprintf(myfile,'bus real appr');
fclose(myfile);

data_voltagem = [ ...
		(1:n-1)'					...
		abs(u_true)				...
		abs(u_appr)];

save('-append', '-ascii', fname, 'data_voltagem');


% Export angles

fname = sprintf('data_%i_voltagea.data', figno);

myfile=fopen(fname,'w');
fprintf(myfile,'bus real appr');
fclose(myfile);

data_voltagea = [ ...
		(1:n-1)' ...
		angle(u_true)/pi*180 ...
		angle(u_appr)/pi*180];

save('-append', '-ascii', fname, 'data_voltagea');


rho_2_star = max(norm(X,2,'cols'));
fprintf(1,'X 2star-norm:   %f\n', rho_2_star);
fprintf(1,'s 2-norm:       %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2_star * norm(s,2));

rho_2 = norm(X,2);
fprintf(1,'X 2-norm:       %f\n', rho_2);
fprintf(1,'s 2-norm:       %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2 * norm(s,2));

rho_1_star = max(max(abs(X)));
fprintf(1,'X 1star-norm:   %f\n', rho_1_star);
fprintf(1,'s 1-norm:       %f\n', norm(s,1));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_1_star * norm(s,1));

rho_Inf_star = norm(X, Inf);
fprintf(1,'X infstar-norm: %f\n', rho_Inf_star);
fprintf(1,'s inf-norm:     %f\n', norm(s,Inf));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_Inf_star * norm(s,Inf));

function lambda = minlambda(y, B)

n = size(B,1);
A = multiprod(y,B,[1],[3]);

%lambda = real(eigs(A,1));
%return

lambda = 0;
step = 0.1;

while step > 1e-5

	[~,p] = chol(lambda * eye(n) - A);

		if p==0

			lambda = lambda - step/2;
			step = step/2;

		else

			lambda = lambda + step;

		end

end

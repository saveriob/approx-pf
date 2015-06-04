lambda = 0;
step = 0.1;

while step > 1e-3

	[~,p] = chol(lambda * eye(n-1) - A);

		if p==0

			lambda = lambda - step/2;
			step = step/2;

		else

			lambda = lambda + step;

		end

end

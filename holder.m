function q = holder(p)

	if p == 1
		q = Inf;
	elseif p == Inf
		q = 1;
	else
		q = p / (p-1);
	end


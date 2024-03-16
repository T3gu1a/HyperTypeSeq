#product of m-fold indicator sequences -- Part II

mfoldCoincidence := proc(L1::[nonnegint, nonnegint], 
			 L2::[nonnegint, nonnegint], $)::[Or(nonnegint, extended_numeric), nonnegint]; 
		local k::nonnegint, mu::posint, J1::set(nonnegint), J2::set(nonnegint), r::nonnegint; 
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "Product of m-fold indicator sequences: finding the coincidence";
		mu := lcm(L1[2], L2[2]); 
		J1 := {seq(L1[1] + k*L1[2], k = 0 .. floor((mu - L1[1])/L1[2]))}; 
		J2 := {seq(L2[1] + k*L2[2], k = 0 .. floor((mu - L2[1])/L2[2]))}; 
		r := min(J1 intersect J2); 
		return [r, mu];
	end proc:
	

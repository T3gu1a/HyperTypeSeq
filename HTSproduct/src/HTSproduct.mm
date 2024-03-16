#Product of hypergeometric-type terms

HTSproduct := proc(s1::algebraic, s2::algebraic, n::name, $)::algebraic; 
		local indetsChis1::list, coefChis1::list, n1::nonnegint, 
		s0::algebraic, s::algebraic, indetsChis2::list, coefChis2::list, 
		n2::nonnegint, i::posint, j::posint, indetsChi::set; 
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "Product of hypergeometric-type sequences";		
		indetsChis1 := [op(indets(s1, 'chi[set(`=`)]'))]; 
		coefChis1 := map(r -> coeff(s1, r), indetsChis1); 
		n1 := numelems(coefChis1); 
		s0 := s1 - add(indetsChis1[i]*coefChis1[i], i = 1 .. n1); 
		s := s0*s2; 
		indetsChis2 := [op(indets(s2, 'chi[set(`=`)]'))]; 
		coefChis2 := map(r -> coeff(s2, r), indetsChis2); 
		n2 := numelems(coefChis2); 
		s0 := s2 - add(indetsChis2[j]*coefChis2[j], j = 1 .. n2); 
		s := s + s0*s1;
		s := s + add(add(coefChis1[i]*coefChis2[j]*mfoldProduct(indetsChis1[i], indetsChis2[j], n), j = 1 .. n2), i = 1 .. n1); 
		indetsChi := indets(s, 'chi[set(`=`)]'); 
		return collect(s, indetsChi, 'distributed'); 
	end proc:
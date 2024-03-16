#product of m-fold indicator sequences -- Part I

mfoldProduct := proc(mf1::Or(chi[set(`=`)], identical(0), identical(1)), 
		     mf2::Or(identical(0), identical(1), chi[set(`=`)]), 
		     n::name, 
		     $)::Or(indentical(0), chi[set(`=`)]); 
			local m1::nonnegint, j1::nonnegint, m2::nonnegint, j2::nonnegint, 
			mu::nonnegint, j::nonnegint; 
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description "Product of m-fold indicator sequences";
			if mf1 = 1 or mf2 = 1 then 
				return mf2*mf1 
			end if; 
			if mf1 = 0 or mf2 = 0 then 
				return 0 
			end if; 
			if n = op([1, 1, 1, 2, 1], mf1) and n = op([1, 1, 1, 2, 1], mf2) then 
				j1 := op([1, 1, 2], mf1); 
				j2 := op([1, 1, 2], mf2); 
				m1 := op([1, 1, 1, 2, 2], mf1); 
				m2 := op([1, 1, 1, 2, 2], mf2); 
				j, mu := op(mfoldCoincidence([j1, m1], [j2, m2])); 
				if j = infinity then 
					return 0 
				else 
					return chi[{modp*(n, mu) = j}]
				end if 
			end if 
		end proc:

#find all representations of m-fold hypergeometric terms for a given m
allhyperterms := proc(hm::list,Cn::list,n::name,algo::procedure,$)
		local m, Lh, j, hypmj;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		m:=hm[1];
		Lh:=[[[op(subs(n=n/m,hm[2]))],[m,0]]];
		for j to m-1 do
			hypmj:=map(allvalues,mfoldterms(Cn,n,[m,j,algo],pochsimpfps));
			Lh:=[op(Lh),[[op(subs(n=(n-j)/m,hypmj))],[m,j]]]
		end do;
		op(Lh)
	end proc:
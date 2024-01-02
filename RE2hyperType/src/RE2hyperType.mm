#two-term recurrence relation case

RE2hyperType := proc(Cn::list,
		n::name,
		s::Or(algebraic,procedure,list),
		chi::name,
		$)::Or(algebraic,identical(FAIL));
		local m::posint, R::ratpoly, hypterms::list, j::nonnegint, V::list, 
		      n0::nonnegint,E0::list(nonnegint),c::nothing, eqs::list, soleqs::list;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "two-term recurrence relation case";
		#There is only one characteristic here 
		#and it is equal to the order of the recurrence
		m:=numelems(Cn)-1;
		#compute the ratio and find the hypergeometric term formulas
		R:=normal(-Cn[1]/Cn[-1]);
		hypterms:=[seq(subs(n=m*n+j,R),j=0..m-1)];
		hypterms:=map(h->pochfactorsimp(h,n,pochsimpfps),hypterms);
		hypterms:=map(h->simplify(h,factorial),[seq(subs(n=(n-j)/m,hypterms[j+1]),j=0..m-1)]);
		#use initial values to determine the linear combination. 
		#V contains the unknowns coefficients
		V:=[seq(c[j],j=0..m-1)];
		n0:=findn0(Cn[1]*Cn[-1],n,m);
		E0:=[seq(j,j=n0..n0+m-1)];
		#linear system encoding the initial conditions
		if type(s,procedure) then
			eqs := [seq(s(E0[j]) = V[j]*eval(hypterms[j],n=E0[j])*mfolder(E0[j],m,j-1),j=1..m)]
		elif type(s,list) then
			eqs := [seq(s[j+1] = V[j+1]*eval(hypterms[j+1],n=j)*mfolder(j,m,j),j=0..m-1)]
		else
			eqs := [seq(try expand(eval(s,n=E0[j])) catch: Limit(s,n=E0[j]) end try 
					 = V[j]*eval(hypterms[j],n=E0[j])*mfolder(E0[j],m,j-1),j=1..m)]
		end if;
		soleqs:=SolveTools:-Linear(eqs, V);
		if soleqs=NULL then return FAIL end if;
		V:=subs(soleqs,V);
		return subs(chi[{`mod`*(n,1)=0}]=1,add(V[j+1]*hypterms[j+1]*chi[{`mod`*(n,m)=j}],j=0..m-1));
	end proc:

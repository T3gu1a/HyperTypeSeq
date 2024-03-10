#November 13, 2023: improving the algorithm for hypergeometric type sequences
#March 10, 2024: creation of HTStoHolonomicRE for finding recurrences from hypergeometric-type normal forms


HolonomicRE := proc(expr::algebraic,
			F::function(name),
			mreord::posint:=10,
			rshift::posint:=1,
			{maxreorder::posint:=10,
			reshift::posint:=1,
			partialwrt::name:=NULL},
			$)::Or(equation,identical(FAIL));
		local partial, Nmax, z, f, N, Coef, Dre, i, mdord, rstp;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description	" Compute a holonomic recurrence equation (mostly of lowest order)"
				" in terms of shifts of F:=F(z), F and z given, satisfied by "
				" a given expression expr depending on z.";
		#small changes
		mdord:=max(maxreorder,mreord);
		rstp:=max(reshift,rshift);
		#multivariate case
		if (partialwrt in {op(F)}) and numelems({op(F)})>1 then
			return PartialHolonomicRE(expr,F,mdord,rstp,case,partialwrt)
		elif numelems({op(F)})>1 then
			return [seq(PartialHolonomicRE(expr,F,mdord,rstp,case,partial),partial in {op(F)})]
		end if;
		#First order case: hyperexponential sequences
		f:=FirstOrderRE(expr,F,rstp);
		if f<>false then
			return f
		end if;
		#Making sure the given maximum DE order is at least equal to 10
		Nmax:=max(10,mdord);
		z:=op(F); 
		#write the ansatz at order 0
		f:=ifelse(has(expr,[cos,cosh,sin,sinh]),
			ansatz(expr,2),ansatz(expr));
		if f=0 then 
			return F=0
		end if;
		#Computing the coefficients (as rational functions) of the RE sought and its order
		Coef,N:=ComputHRE(f,z,Nmax,rstp);
		#if Coef is not empty then build the DE
		if Coef<>[] then
			#clearing denominators
			Dre:=lcm(op(denom(Coef)));
			#multiply all coefficients by the lcm of their denominators
			Coef:=map(r-> factor(normal(Dre*r)), Coef);
			Coef:=[op(Coef),Dre];
			#The zeroth order is handled with the empty list (see the diff help page)
			return add(Coef[i+1]*LREtools:-shift(F,z,i*rstp),i=0..N)=0	
		else
			#No holonomic DE of order less than Nmax+1 is found
			userinfo(2,HolonomicRE,printf("No holonomic RE of order less than %d found\n", Nmax+1));
			return FAIL
		end if
	end proc:

HTStoHolonomicRE:= proc(expr::algebraic,a::anyfunc(name),{addorder::nonnegint:=0,restep::posint:=1})::equation;
		local indetsChi::list,coefChi::list,n::name:=op(a),REs::list,chi0::algebraic,
		      r::nothing,u::nothing,i::posint,j::nonnegint;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "Compute the recurrence of hypergeometric-type term";
		indetsChi:=[op(select(has,indets(expr),chi))];
		coefChi:=map(r-> [coeff(expr,r),r], indetsChi);
		chi0:=expr-add(r[1]*r[2],r in coefChi);
		if chi0 <> 0 then
			coefChi:=[op(coefChi),[chi0,chi[modp*(n,1)=0]]];
		end if;
		REs:=map(r->[r[1],op([1,1,2,2],r[2]),op([1,2],r[2])],coefChi);
		REs:=map(r->[subs(n=r[2]*n+r[3],r[1]),r[2],r[3]],REs);
		REs:=map(r->[HolonomicRE(r[1],u(n)),r[2],r[3]],REs);
		REs:=map(r->[subs([seq(u(n+j)=u(r[2]*n+r[3]+r[2]*j),j=0..REorder(r[1],u(n)))],r[1]),r[2],r[3]],REs);
		REs:=map(r-> normal(subs(n=(n-r[3])/r[2],r[1])),REs);
		REs:=map(r->REcoeff(r,u(n)),REs);
		REs:=map(simpRE,REs);
		REs:=[seq(add(REs[i][j+1]*u[i](n+j),j=0..numelems(REs[i])-1)=0,i=1..numelems(REs))];
		AddHolonomicRE(REs,[seq(u[i](n),i=1..numelems(REs))],a,':-addorder'=addorder,':-restep'=restep)
	end proc:
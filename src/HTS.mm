#The exported algorithms

#convert analogue
`convert/HypergeometricType` := proc(s::Or(algebraic,`=`),
				n::Or(name,[anyfunc(name),procedure]),
				{maxreorder::posint:=10,
				varreq::Or(name,identical(NULL)):=NULL,
				mfoldvar::Or(name,identical(NULL)):=NULL},
				$)::Or(algebraic,identical(FAIL));
		option `Copyright (c) 2023 Bertrand Teguia T.`;			
		if type(s,`=`) then
			return REtoHTS(s,n[1],n[2],':-mfoldvar'=mfoldvar)
		else
			return HTS(s,n,':-maxreorder'=maxreorder,':-varreq'=varreq,':-mfoldvar'=mfoldvar)
		end if
	end proc:
	
REtoHTS := proc(RE::`=`,
		S::anyfunc(name),
		F::Or(procedure,list),
		{mfoldvar::Or(name,identical(NULL)):=NULL},
		$)::Or(algebraic,identical(FAIL));
		local Cn::list, n::name, chi::name, hypertype::algebraic, d::posint, j::nonnegint;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "Closed form of hypergeometric type sequences from a given recurrence relation";
		#interlacement variable
		if mfoldvar=NULL then 
			`tools/genglobal`('chi',{},'reset');
			chi:=`tools/genglobal`('chi')
		else
			chi:=mfoldvar
		end if;
		n:=op(S);
		#collect the coefficients
		Cn:=REcoeff(lhs(RE)-rhs(RE)=0,S);
		d:=numelems(Cn)-1;
		if type(F,list) then
			userinfo(2,HTS,printf("Not enough initial values\n"));
			if numelems(F)<d then
				return FAIL
			end if
		end if;
		#the recurrence equation has only two nonzero terms
		if numelems(select(nonzerop,Cn))=2 then
			return RE2hyperType(Cn,n,F,chi)
		end if;
		#Use the general hypergeometric type algorithm
		hypertype:=HyperTypeForm(Cn,n,F,d,chi);
		if hypertype<>false then
			return hypertype
		else 
			#No hypergeometric type formula found
			userinfo(2,HTS,printf("No hypergeometric type formula found\n"));
			if type(F,procedure) then
				return RE,seq(eval(S,n=j)=F(j),j=0..d-1)
			else
				return RE,seq(eval(S,n=j)=F[j+1],j=0..d-1)
			end if
		end if;
	end proc:

#main procedure
HTS := proc(s::algebraic,
		n::name,
		{maxreorder::posint:=10,
		varreq::Or(name,identical(NULL)):=NULL,
		mfoldvar::Or(name,identical(NULL)):=NULL},
		$)::Or(algebraic,identical(FAIL));
		local S::name, chi::name, hypertype::algebraic, RE::algebraic, Cn::list, d::posint, j::nonnegint;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "Closed form of hypergeometric type sequences from a given expression";
		#compute the recurrence equation
		#the global variable S will appear in the output (recursive formula) if
		#no hypergeometric type series is found
		if varreq=NULL then 
			`tools/genglobal`('S',{},'reset');
			S:=`tools/genglobal`('S')
		else
			S:=varreq
		end if;
		#interlacement variable
		if mfoldvar=NULL then 
			`tools/genglobal`('chi',{},'reset');
			chi:=`tools/genglobal`('chi')
		else
			chi:=mfoldvar
		end if;
		RE:=HolonomicRE(s,S(n),maxreorder);
		if RE<>FAIL then
			#collect the coefficients
			Cn:=REcoeff(RE,S(n));
			#the recurrence equation has only two nonzero terms
			if numelems(select(nonzerop,Cn))=2 then
				return RE2hyperType(Cn,n,s,chi)
			end if;
			#Use the general hypergeometric type algorithm
			hypertype:=HyperTypeForm(Cn,n,s,maxreorder,chi);
			if hypertype<>false then
				return hypertype
			else
				#No hypergeometric type formula found
				userinfo(2,HTS,printf("No hypergeometric type formula found\n"));
				d:=numelems(Cn);
				return RE,seq(S(j)=eval(s,n=j),j=0..d-1)
			end if;
		else
			userinfo(2,HTS,printf("No holonomic RE of order at most %d found\n", maxreorder));
			return FAIL
		end if
	end proc:

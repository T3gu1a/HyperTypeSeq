#General hypergeometric type case
HyperTypeForm := proc(Cn::list,
		n::name,
		s::Or(algebraic,procedure,list),
		REorder::posint,
		chi::name,
		$)::Or(algebraic,indentical(false));
		local algo::procedure, Hyp::list(algebraic), A::nothing, 
		      CCn::list, t::posint, hyperType::algebraic;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		algo:=HyperTermQ;
		#compute m-fold hypergeometric term solutions. Cn is the list of coefficients
		Hyp:=try mfoldterms(Cn,n,[algo],pochsimpfps) catch: [] end try;
		CCn:=copy(Cn);
		#if there are no solutions
		if Hyp=[] then
			#Try computations from 2-fold and 3-fold holonomic REs
			#this is to avoid representations with complex numbers as much as possible
			if numelems(indets(Cn,'name'))<3 and type(s,algebraic) then 
				t:=2;
				while Hyp=[] and t<4 do
					CCn:=REcoeff(HolonomicRE(s,A(n),REorder+t,t),A(n));
					if numelems(select(nonzerop,CCn))=2 then
						return RE2hyperType(CCn,n,s,chi)
					end if;
					Hyp:=try mfoldterms(CCn,n,[algo],pochsimpfps) catch: [] end try;
					t:=t+1
				end do;
			end if;
			#If there are no solution then consider the solutions
			# over algebraic extension fields by computing with HyperTermC
			if Hyp=[] then
				algo:=HyperTermC;
				CCn:=copy(Cn);
				Hyp:=map(r->[r[1],map(allvalues,r[2])],mfoldterms(Cn,n,[algo],pochsimpfps));
			end if
		end if;
		if Hyp<>[] then
			#Find the linear combination of hypergeometric type sequences
			hyperType:=HyperTypeCombin(Hyp,algo,CCn,n,s,chi);
			if hyperType<>false then
				return hyperType
			else
				#otherwise if the algebraic case was not applied then we apply it
				if algo=HyperTermQ then
					algo:=HyperTermC;
					Hyp:=map(r->[r[1],map(allvalues,r[2])],mfoldterms(CCn,n,[algo],pochsimpfps));
					hyperType:=HyperTypeCombin(Hyp,algo,CCn,n,s,chi);
					return ifelse(hyperType=false,false,hyperType)
				end if
			end if
		end if;
		return false
	end proc:

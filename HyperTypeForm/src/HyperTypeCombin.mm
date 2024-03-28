#Find the linear combination of hypergeometric type sequences
HyperTypeCombin := proc(Hyp::list, 
		  algo::procedure,
		  Cn::list,
		  n::name,
		  s::Or(algebraic,procedure,list),
		  chi::name,
		  $)::Or(algebraic,identical(false));
		local  mftypes::list, M::posint, N::posint, V::list, Vc::list, c::nothing, soleqs::list,
		      n0::nonnegint,E0::list(nonnegint),i::posint,j::nonnegint, k::posint, mft, eqs::list(algebraic);
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		#Complete list of all m-fold hypergeometric terms for each m leading to solutions.
		#Here all the necessary data for the hypergeometric type representations are provided.
		mftypes:=map(hm->allhyperterms(hm,Cn,n,algo),Hyp);
		M:=numelems(mftypes);
		N:=add(numelems(mftypes[j][1]),j=1..M);
		N:=max(N,numelems(Cn));
		V:=[seq(c[j],j=1..N)];
		j:=0;
		i:=1;
		Vc:=copy(V);
		#Evaluations of hypergeometric type series for finding the linear combination
		for mft in mftypes do
			mftypes[i]:=[[add(V[j+k]*convert(simplify(mft[1][k],factorial),factorial),k=1..numelems(mft[1]))],mft[2]];
			j:=j+numelems(mft[1]);
			i:=i+1
		end do;
		n0:=findn0(Cn[1]*Cn[-1],n,N);
		E0:=[seq(j,j=n0..n0+N-1)];
		if type(s,procedure) then
			eqs:=[seq(s(E0[j])=add(eval(op(mftypes[k][1]),n=E0[j])*mfolder(E0[j],op(mftypes[k][2])),k=1..M),j=1..N)]
		elif type(s,list) then	
			eqs:=[seq(s[j+1]=add(eval(op(mftypes[k][1]),n=j)*mfolder(j,op(mftypes[k][2])),k=1..M),j=0..N-1)]
		else
			eqs:=[seq(try expand(HTSeval(s,n=E0[j])) catch: Limit(s,n=E0[j]) end try
				=add(eval(op(mftypes[k][1]),n=E0[j])*mfolder(E0[j],op(mftypes[k][2])),k=1..M),j=1..N)]
		end if;
		soleqs:=SolveTools:-Linear(eqs, V);
		Vc:=map(v->v=0,Vc);
		V:=subs(soleqs,V);
		V:=subs(Vc,V);
		V:=seq(lhs(Vc[j])=V[j],j=1..N);
		#writing of the linear combination
		if soleqs<>NULL then
			mftypes:=subs(V,mftypes);
			return subs(chi[{`mod`*(n,1)=0}]=1,add(op(mftypes[k][1])*chi[{`mod`*(n,mftypes[k][2][1])=mftypes[k][2][2]}],k=1..M))
		end if;
		return false
	end proc:

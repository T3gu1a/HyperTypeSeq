#adding, multiplying, and exponentiating holonomic sequences from their recurrence equations
	
FindCoefsHRE := proc(f::algebraic,
			z::name,
			Nmax::posint,
			restep::posint,
			Subshift::list(`=`),
			startfromord::posint:=1,
			$)
		local df, N:=0, term, Basis, Comp, Mat, Nb, Coef:=[];
		option `Copyright (c) 2020 Bertrand Teguia T.`;
		df:=f;
		if df=0 then
			return [1], N
		end if;
		term:=colTerm(df,`+`);
		Basis,Comp:=augBasis(term,z);
		Mat:=Matrix(Comp);
		Nb:=numelems(Basis);
		while Coef=[] and N<Nmax do
			N:=N+1; 
			df:=LREtools:-shift(df,z,restep);
			to restep do
				df:=eval(df,Subshift)
			end do;
			df:=expand(df);
			term:=colTerm(df,`+`);
			Basis,Comp:=basisComp(Basis,term,z);
			Mat:=adapt(Mat,Comp);
			if Nb=numelems(Comp) then
				Coef:=holoSolve(Mat,N)
			else
				Nb:=numelems(Comp)
			end if
		end do;
		return Coef, N
	end proc:

SystemtoHRE := proc(f::algebraic,
			F::anyfunc(name),
			sublistshift::list(`=`),
			{maxreorder::posint:=10,
			restep::posint:=1,
			startfromord::posint:=1},
			$)::Or(equation,identical(FAIL));
		local Nmax, z, N, Coef, Dde, i;
		option `Copyright (c) 2020 Bertrand Teguia T.`;
		description " to be...";
		Nmax:=max(10,startfromord+maxreorder);
		z:=op(1,F);
		Coef,N:=FindCoefsHRE(f,z,Nmax,restep,sublistshift,startfromord);
		if Coef<>[] then
			Dde:=lcm(op(denom(Coef)));
			Coef:=map(r-> factor(normal(Dde*r)), Coef);
			Coef:=[op(Coef),factor(Dde)];
			return add(Coef[i+1]*LREtools:-shift(F,z,i*restep),i=0..N)=0	
		else
			userinfo(2,HolonomicRE,printf("No holonomic DE of order less than %d found\n", Nmax+1));
			return FAIL
		end if
end proc:

buildHRESystem:= proc(DE::`=`,
			y::anyfunc(name),
			x::name,
			$)::list(`=`);
	local d::posint, t::name, SubL::list, PolDE::polynom, j::posint;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	t:=op(y);
	d:=REorder(DE,y);
	SubL:=[seq(LREtools:-shift(y,t,j)=x[j],j=0..d)];
	PolDE:=subs(SubL,lhs(DE));
	[[seq(x[j],j=1..(d-1)),solve(PolDE,x[d])],[seq(x[j],j=0..(d-1))]]
end proc:

mergeHRESystem:= proc(L::list(`=`),
			V::list(anyfunc(name)),
			$)::`=`;
	local l::posint:=numelems(L), j::posint, Sys::list, vars::list, deriv::list, 
	      n::posint, x::nothing, X::list, i::posint, Ind::list;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	Sys:=[seq(buildHRESystem(L[j],V[j],cat(x,j)),j=1..l)];
	vars:=map(r->op(r[2]),Sys);
	deriv:=map(r->op(r[1]),Sys);
	n:=numelems(vars);
	X:=[seq(vars[j]=x[j],j=1..n)];
	Ind:=[seq(1+add(numelems(Sys[i][2]),i=1..(j-1)),j=1..l)];
	[subs(X,deriv),map(r->x[r],Ind),map(rhs,X)]
end proc:

AddHolonomicRE:= proc(L::list(`=`),
			   V::list(anyfunc(name)),
			   z::anyfunc(name),
			   {addorder::nonnegint:=0,
			   restep::posint:=1},
			   $)::Or(`=`,FAIL);
	local t::name, start::posint, reords::list, endord::posint, j::posint, Sys::list, subvars::list, SubL::list;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	t:=op(z);
	reords:=[seq(REorder(L[j],V[j]),j=1..numelems(L))];
	start:=min(reords); 
	endord:=add(reords);
	if numelems(L)=1 then
		return subs(op(0,V[1])=op(0,z),L[1])
	end if;
	Sys:=mergeHRESystem(L,V);
	subvars:=map(r->r=r(t),Sys[3]);
	Sys:=subs(subvars,Sys);
	SubL:=[seq(LREtools:-shift(Sys[3][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
	SystemtoHRE(add(Sys[2]),z,SubL,maxreorder=endord+addorder,startfromord=start,':-restep'=restep)
end proc:

MulHolonomicRE:= proc(L::list(`=`),
			   V::list(anyfunc(name)),
			   z::anyfunc(name),
			   {addorder::nonnegint:=0,
			   restep::posint:=1},
			   $)::Or(`=`,FAIL);
	local t::name, start::posint, reords::list, endord::posint, j::posint, Sys::list, subvars::list, SubL::list;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	t:=op(z);
	reords:=[seq(REorder(L[j],V[j]),j=1..numelems(L))];
	start:=min(reords); 
	endord:=add(reords);
	if numelems(L)=1 then
		return subs(op(0,V[1])=op(0,z),L[1])
	end if;
	Sys:=mergeHRESystem(L,V);
	subvars:=map(r->r=r(t),Sys[3]);
	Sys:=subs(subvars,Sys);
	SubL:=[seq(LREtools:-shift(Sys[3][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
	SystemtoHRE(mul(Sys[2]),z,SubL,maxreorder=endord+addorder,startfromord=start,':-restep'=restep)
end proc:

SelfOpHolonomicRE:= proc(DE::`=`,
			y::anyfunc(name),
			z::name=polynom,
			{addorder::posint:=10,
			restep::posint:=1},
			$)::Or(`=`,FAIL);
	local t::name:=op(y), start::posint, Sys::list,x::nothing,var::name,subvars::list,SubL::list,j::posint;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	start:=REorder(DE,y);
	Sys:=buildHRESystem(DE,y,x);
	var:=op(0,y);
	subvars:=map(r->r=r(t),Sys[2]);
	Sys:=subs(subvars,Sys);
	SubL:=[seq(LREtools:-shift(Sys[2][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
	SystemtoHRE(subs(var=Sys[2][1],normal(rhs(z))),lhs(z)(t),SubL,maxreorder=start+addorder,startfromord=start,':-restep'=restep)
end proc:

REorder := proc(RE::Or(algebraic,`=`),a::anyfunc(name),$)::nonnegint;
		local n::name,aterms::list,A;
		option `Copyright (c) 2020 Bertrand Teguia T.`;
		n:=op(1,a);
		A:=op(0,a);
		aterms:=indets(RE,A(`+`)) union indets(RE,A(name));
		aterms:=subs(n=0,map(r-> op(r), aterms));
		return max(aterms)
	end proc:
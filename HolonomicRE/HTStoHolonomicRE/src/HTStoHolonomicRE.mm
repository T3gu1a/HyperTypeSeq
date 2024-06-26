

#March 10, 2024: creation of HTStoHolonomicRE for finding recurrences from hypergeometric-type normal forms

HTStoHolonomicRE:= proc(expr::algebraic,a::anyfunc(name),{addorder::nonnegint:=0,restep::posint:=1})::`=`;
		local indetsChi::list,coefChi::list,n::name:=op(a),REs::list,chi0::algebraic,
		      r::nothing,u::nothing,i::posint,j::nonnegint;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "Compute the recurrence of hypergeometric-type term";
		indetsChi:=[op(select(has,indets(expr),chi))];
		coefChi:=map(r-> [coeff(expr,r),r], indetsChi);
		chi0:=expr-add(r[1]*r[2],r in coefChi);
		if chi0 <> 0 then
			coefChi:=[op(coefChi),[chi0,chi[{modp*(n,1)=0}]]];
		end if;
		REs:=map(r->[r[1],op([1,1,1,2,2],r[2]),op([1,1,2],r[2])],coefChi);
		REs:=map(r->[subs(n=r[2]*n+r[3],r[1]),r[2],r[3]],REs);
		REs:=map(r->[HolonomicRE(r[1],u(n)),r[2],r[3]],REs);
		REs:=map(r->[subs([seq(u(n+j)=u(r[2]*n+r[3]+r[2]*j),j=0..REorder(r[1],u(n)))],r[1]),r[2],r[3]],REs);
		REs:=map(r-> normal(subs(n=(n-r[3])/r[2],r[1])),REs);
		REs:=map(r->REcoeff(r,u(n)),REs);
		REs:=map(simpRE,REs);
		REs:=[seq(add(REs[i][j+1]*u[i](n+j),j=0..numelems(REs[i])-1)=0,i=1..numelems(REs))];
		AddHolonomicRE(REs,[seq(u[i](n),i=1..numelems(REs))],a,':-addorder'=addorder,':-restep'=restep)
	end proc:
	

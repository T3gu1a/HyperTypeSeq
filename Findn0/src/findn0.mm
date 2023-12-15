#Finding the set of integers to evaluate
#the ansatz of the linear combination

findn0:=proc(P::algebraic,n::name,N::posint)::nonegint;
	local n0::nonnegint:=0,j::nonnegint,L::list;
	L:=[seq(eval(P,n=j),j=0..N)];
	while 0 in L do:
		n0:=n0+[ListTools:-SearchAll(0,L)][-1];
		L:=[seq(eval(P,n=j),j=n0..n0+N)]
	end do;
	return n0
end proc:

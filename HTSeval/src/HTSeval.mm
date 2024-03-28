HTSeval := proc(s::algebraic,eq::name=nonnegint)::algebraic;
		local indetsChis::list, i::postint, mjs::list, 
		n::name:=lhs(eq), k::nonnegint:=rhs(eq), Sub::list;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		indetsChis := [op(indets(s, 'chi[set(`=`)]'))];
		mjs := map(r-> [op([1, 1, 1, 2, 2], r), op([1, 1, 2], r)],indetsChis);
		Sub := [seq(indetsChis[i]=mfoldInd(k,mjs[i][1],mjs[i][2]),i=1..numelems(mjs))];
		return eval(subs(Sub,s),n=k) 
	end proc:

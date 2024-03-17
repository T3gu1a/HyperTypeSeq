#procedures for interlacements when evaluating the sequence
mfolder := proc(n::nonnegint, m::posint, j::nonnegint,$)::Or(indentical(0),identical(1)); 
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "procedure for valued interlacements during evaluation";
		if n mod m = j mod m then 
			return 1
		else 
			return 0
		end if
	end proc:

#procedure for writing m-fold indicator terms
mfoldInd := proc(n::polynom,m::nonnegint,j::nonnegint,$)
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "the m-fold indicator term";
		if m=0 then 
			return 0
		end if;
		if m=1 then
			return 1
		end if;
		if type(n mod m, nonnegint) then
			mfolder(n,m,j)
		else
			return chi[{modp*(n,m) = j}]
		end if
	end proc:
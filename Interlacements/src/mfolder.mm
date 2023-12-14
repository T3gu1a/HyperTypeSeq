#procedure for interlacements
mfolder := proc(n::nonnegint, m::posint, j::nonnegint,$)::Or(indentical(0),identical(1)); 
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "two-term recurrence relation case";
		if n mod m = j mod m then 
			return 1
		else 
			return 0
		end if
	end proc:
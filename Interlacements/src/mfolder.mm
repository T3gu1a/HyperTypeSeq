#procedure for interlacements
mfolder := proc(n::nonnegint, m::posint, j::nonnegint,$)::Or(indentical(0),identical(1)); 
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "m-fold indicator sequence";
		if n mod m = j mod m then 
			return 1
		else 
			return 0
		end if
	end proc:
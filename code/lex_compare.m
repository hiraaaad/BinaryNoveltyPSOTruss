function result = lex_compare(x,y)
    result = false;
   % true means that x is better that y
   if x(1) <= y(1)
       result = true;
       return
   elseif x(1) == y(1) && x(2) <= y(2)
        result = true;
        return
   elseif x(1) == y(1) && x(2) == y(2) && x(3) <= y(3)
        result = true;
        return
   else
       result = false;
   end
end
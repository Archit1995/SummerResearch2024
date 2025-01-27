
a:=36;
f:=ifactor(a);

nops(op(1,(op(2,f))));

prime_var_map:=table([3=x,5=y,7=z,11=w]);

m1:=subs(op(1,(op(1,f)))=x,f);
m2:=subs(op(1,(op(2,m1)))=y,m1);

c:=1155*3*25;
monomial_generator:=proc(roots_,prime_var_map)
     local ff,l,l2,i;
     monomials:=Vector(numelems(roots_,0));
     for j in numelems(roots_) do 
          ff:=ifactor(roots_[j]);
          l:=nops(ff);
          for i from 1 to l do
               l2:=nops(op(i,ff)):
               if l2=1 then 
                    ff:=subs(op(i,ff)=prime_var_map[op(1,op(i,ff))],ff);
               else 
                    ff:=subs(op(1,op(i,ff))=prime_var_map[op(1,(op(1,op(i,ff))))],ff);
               fi;
          end do;
          monomials[j]:=ff;
     end do;
     return convert(monomials,list);
end proc;
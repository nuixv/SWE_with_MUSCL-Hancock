% minmod function
       function m=minmod(a,b)
          if (abs(a)<=abs(b))&&(a*b>0)
             m=a;
          elseif (abs(b)<=abs(a))&&(a*b>0)
             m=b;
          elseif (a*b<=0)
             m=0;
          end
       end

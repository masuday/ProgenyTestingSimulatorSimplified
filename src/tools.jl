function extract_data(y)
   nrec = length(y)
   if nrec>0
      ntr = length(y[1])
   else
      ntr = 0
   end
   data = zeros(ntr,nrec)
   for j=1:nrec
      for i=1:ntr
         data[i,j] = y[j][i]
      end
   end
   return data
end

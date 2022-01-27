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

"""
    fix_covariance_matrix!(M, value=1e-6)

Fix covariance matrix `M` to be positive definite by replacing negative eigenvalues with a small number.
"""
function fix_covariance_matrix!(M, value=1e-6)
   EVEC = eigen(M).vectors
   EVAL = diagm(eigen(M).values)
   EVAL[EVAL .< 0.0] .= value
   M .= EVEC * EVAL * EVEC'
end

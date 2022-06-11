"""
    write_data_1st(io::IO, df::DataFrame)

Write a data file for BLUPF90 programs assuming a simple animal model with the 1st lactation records.
"""
function write_data_1st(io::IO, df::DataFrame)
   nrec = size(df,1)
   for i=1:nrec
      if !ismissing(df.y[i][1])
         print(io, @sprintf("%d %d %.2f\n",i,df.herd[i],df.y[i][1]))
      end
   end
end

"""
    write_data_rep(io::IO, df::DataFrame)

Write a data file for BLUPF90 programs assuming a repeatability model with all the records available.
"""
function write_data_rep(io::IO, df::DataFrame)
   nrec = size(df,1)
   nlac = length(df.y[1])
   for i=1:nrec
      for k=1:nlac
         if !ismissing(df.y[i][k])
            print(io, @sprintf("%d %d %.2f\n",i,df.herd[i],df.y[i][k]))
         end
      end
   end
end

"""
    write_data_mt(io::IO, df::DataFrame; ntraits)

Write a data file for BLUPF90 programs assuming a multiple-trait model with all the records available.
"""
function write_data_mt(io::IO, df::DataFrame; ntraits::Union{Nothing,Int64}=nothing)
   nrec = size(df,1)
   if isnothing(ntraits)
      nlac = length(df.y[1])
   else
      nlac = max(1,min(ntraits,length(df.y[1])))
   end
   for i=1:nrec
      nmiss = sum(ismissing.(df.y[i]))
      if nmiss < nlac
         print(io, @sprintf("%d %d",i,df.herd[i]))
         for k=1:nlac
            y = ifelse(ismissing(df.y[i][k]), 0.0, df.y[i][k])
            print(io, @sprintf(" %.2f",y))
         end
         print(io, @sprintf("\n"))
      end
   end
end

"""
    write_pedigree(io::IO, df::DataFrame)

Write a pedigree file for BLUPF90 programs with all the animals available in the pedigree.
The following items will be written.

1. animal ID
2. sire ID
3. dam ID
4. inbcode
5. year og birth
6. male flag
7. proven flag
8. number of records
9. number of daughters with records
10. inbreeding coefficient
11. herd
"""
function write_pedigree(io::IO, df::DataFrame)
   nped = size(df,1)
   for i=1:nped
      aid = df.aid[i]
      sid = df.sid[i]
      did = df.did[i]
      inbcode = get_inbupg_code(sid,did,df.f)
      f = df.f[aid]
      gen = df.gen[aid]
      nrec = sum( .!ismissing.(df.y[aid]) )
      nrecdau = df.nrecdau[aid]
      male = df.male[aid] + 0
      proven = df.proven[aid] + 0
      herd = df.herd[aid]
      print(io, @sprintf("%d %d %d %4d %d  %d %d  %d %d %8.5f  %d\n",aid,sid,did,inbcode,gen, male,proven, nrec,nrecdau, f, herd))
   end
end

function get_inbupg_code(s::Int,d::Int,inb::Vector{Float64})
   if s>0
      ms = 0.0
      fs = inb[s]
   else
      ms = 1.0
      fs = 0.0
   end
   if d>0
      md = 0.0
      fd = inb[d]
   else
      md = 1.0
      fd = 0.0
   end
   return round( 4000/((1+ms)*(1-fs) + (1+md)*(1-fd)) )
end

"""
    write_parfile_1st(io::IO, df::DataFrame)

Write a parameter file for BLUPF90 programs assuming a simple animal model with the 1st lactation records.
"""
function write_parfile_1st(io::IO, df::DataFrame, datafile::String, pedfile::String, vg::Float64, ve::Float64; option=Vector{String}[])
   nherd = maximum(df.herd)
   print(io,"# parameter file for a simple animal model
DATAFILE
  $(datafile)
NUMBER_OF_TRAITS
  1
NUMBER_OF_EFFECTS
  2
OBSERVATION(S)
  3
WEIGHT(S)
 
EFFECTS:
  1 $(size(df,1)) cross
  2 $(nherd) cross
RANDOM_RESIDUAL VALUES
  $(ve)
RANDOM_GROUP
  1
RANDOM_TYPE
  add_an_upginb
FILE
  $(pedfile)
(CO)VARIANCES
  $(vg)
#OPTION use_yams
#OPTION EM-REML 5
")
   if length(option)>0
      for str in option
         print(io,str,"\n")
      end
   end
end

"""
    write_parfile_rep(io::IO, df::DataFrame)

Write a parameter file for BLUPF90 programs assuming a repeatability model with all the records available.
"""
function write_parfile_rep(io::IO, df::DataFrame, datafile::String, pedfile::String, vg::Float64, vp::Float64, ve::Float64; option=Vector{String}[])
   nherd = maximum(df.herd)
   print(io,"# parameter file for repeatability model
DATAFILE
  $(datafile)
NUMBER_OF_TRAITS
  1
NUMBER_OF_EFFECTS
  3
OBSERVATION(S)
  3
WEIGHT(S)
 
EFFECTS:
  1 $(size(df,1)) cross
  1 $(size(df,1)) cross
  2 $(nherd) cross
RANDOM_RESIDUAL VALUES
  $(ve)
RANDOM_GROUP
  1
RANDOM_TYPE
  add_an_upginb
FILE
  $(pedfile)
(CO)VARIANCES
  $(vg)
RANDOM_GROUP
  2
RANDOM_TYPE
  diagonal
FILE
 
(CO)VARIANCES
  $(vp)
#OPTION use_yams
#OPTION EM-REML 5
")
   if length(option)>0
      for str in option
         print(io,str,"\n")
      end
   end
end

"""
    write_parfile_mt(io::IO, df::DataFrame, datafile, pedfile, G, R; ntraits=size(G,1))

Write a parfile file for BLUPF90 programs assuming a repeatability model with all the records available.
"""
function write_parfile_mt(io::IO, df::DataFrame, datafile::String, pedfile::String, G::Matrix{Float64}, R::Matrix{Float64}; option=Vector{String}[], ntraits=size(G,1))
   nherd = maximum(df.herd)
   ntr = min(max(1,ntraits),size(G,1))
   effect_line1 = repeat("1 ",ntr) * string(size(df,1)) * " cross"
   effect_line2 = repeat("2 ",ntr) * string(nherd) * " cross"
   obs_line = join(string.([i for i=3:(ntr+2)]), " ")

   print(io,"# parameter file for multiple-trait model
DATAFILE
  $(datafile)
NUMBER_OF_TRAITS
  $(ntr)
NUMBER_OF_EFFECTS
  2
OBSERVATION(S)
  $(obs_line)
WEIGHT(S)
 
EFFECTS:
  $(effect_line1)
  $(effect_line2)
RANDOM_RESIDUAL VALUES
")
   write_matrix(io,R[1:ntr,1:ntr])
   print(io,"RANDOM_GROUP
  1
RANDOM_TYPE
  add_an_upginb
FILE
  $(pedfile)
(CO)VARIANCES
")
   write_matrix(io,G[1:ntr,1:ntr])
#   print(io,"OPTION use_yams
#OPTION EM-REML 10
#")

   if length(option)>0
      for str in option
         print(io,str,"\n")
      end
   end
end

function write_matrix(io::IO, M::Matrix{Float64})
   nr = size(M,1)
   nc = size(M,2)
   for i=1:nr
      for j=1:nc
         #print(io,@sprintf(" %.5f",M[i,j]))
         print(io,@sprintf(" %.8e",M[i,j]))
      end
      print(io,@sprintf("\n"))
   end
end

"""
    load_solutions_1st!(solfile, ebv::Vector{Float64})
    load_solutions_1st!(solfile, ebv::Vector{Float64}, sep::Vector{Float64})

Read EBV from the solution file `solfile` generated by BLUPF90 under a sinple animal model.
The vector `ebv` (with the length of the number of animals) should be prepared by user.
The standard error of prediction, i.e., sqrt(PEV), can be obtained by `sep` if the solution file has the fifth field.
"""
function load_solutions_1st!(solfile, ebv::Vector{Float64})
   load_solutions_rep!(solfile, ebv)
end
function load_solutions_1st!(solfile, ebv::Vector{Float64}, sep::Vector{Float64})
   load_solutions_rep!(solfile, ebv, sep)
end

"""
    load_solutions_rep!(solfile, ebv::Vector{Float64})
    load_solutions_rep!(solfile, ebv::Vector{Float64}, sep::Vector{Float64})

Read EBV from the solution file `solfile` generated by BLUPF90 under a repeatabiliyt model.
The vector `ebv` (with the length of the number of animals) should be prepared by user.
The standard error of prediction, i.e., sqrt(PEV), can be obtained by `sep` if the solution file has the fifth field.
"""
function load_solutions_rep!(solfile, ebv::Vector{Float64})
   open(solfile,"r") do ios
      readline(ios)
      while(!eof(ios))
         line = strip(readline(ios))
         col = split(line)
         trt = parse(Int64,col[1])
         eff = parse(Int64,col[2])
         lev = parse(Int64,col[3])
         sol = parse(Float64,col[4])
         if eff!=1
            break
         end
         if trt>1
            throw(ErrorException("seems not to be repeatability model"))
         end
         if lev>length(ebv)
            throw(DimensionMismatch("short in elements of ebv"))
         else
            ebv[lev] = sol
         end
      end
   end
end

function load_solutions_rep!(solfile, ebv::Vector{Float64}, sep::Vector{Float64})
   open(solfile,"r") do ios
      readline(ios)
      while(!eof(ios))
         line = strip(readline(ios))
         col = split(line)
         trt = parse(Int64,col[1])
         eff = parse(Int64,col[2])
         lev = parse(Int64,col[3])
         sol = parse(Float64,col[4])
         err = parse(Float64,col[5])
         if eff!=1
            break
         end
         if trt>1
            throw(ErrorException("seems not to be repeatability model"))
         end
         if lev>length(ebv)
            throw(DimensionMismatch("short in elements of ebv"))
         else
            ebv[lev] = sol
         end
         if lev>length(sep)
            throw(DimensionMismatch("short in elements of sep"))
         else
            sep[lev] = err
         end
      end
   end
end

"""
    load_solutions_mt!(solfile, ebv::Matrix{Float64})
    load_solutions_mt!(solfile, ebv::Matrix{Float64}, sep::Matrix{Float64})

Read EBV from the solution file `solfile` generated by BLUPF90 under a multiple-trait model.
The matrix `ebv` (the number of traits by the number of animals) should be prepared by user.
The standard error of prediction, i.e., sqrt(PEV), can be obtained by `sep` if the solution file has the fifth field.
"""
function load_solutions_mt!(solfile, ebv::Matrix{Float64})
   open(solfile,"r") do ios
      readline(ios)
      while(!eof(ios))
         line = strip(readline(ios))
         col = split(line)
         trt = parse(Int64,col[1])
         eff = parse(Int64,col[2])
         lev = parse(Int64,col[3])
         sol = parse(Float64,col[4])
         if eff!=1
            break
         end
         if lev>size(ebv,2)
            throw(DimensionMismatch("short in elements for levels of ebv"))
         elseif trt>size(ebv,1)
            throw(DimensionMismatch("short in elements for traits of ebv"))
         else
            ebv[trt,lev] = sol
         end
      end
   end
end

function load_solutions_mt!(solfile, ebv::Matrix{Float64}, sep::Matrix{Float64})
   open(solfile,"r") do ios
      readline(ios)
      while(!eof(ios))
         line = strip(readline(ios))
         col = split(line)
         trt = parse(Int64,col[1])
         eff = parse(Int64,col[2])
         lev = parse(Int64,col[3])
         sol = parse(Float64,col[4])
         err = parse(Float64,col[5])
         if eff!=1
            break
         end
         if lev>size(ebv,2)
            throw(DimensionMismatch("short in elements for levels of ebv"))
         elseif trt>size(ebv,1)
            throw(DimensionMismatch("short in elements for traits of ebv"))
         else
            ebv[trt,lev] = sol
         end
         if lev>size(sep,2)
            throw(DimensionMismatch("short in elements for levels of sep"))
         elseif trt>size(sep,1)
            throw(DimensionMismatch("short in elements for traits of sep"))
         else
            sep[trt,lev] = err
         end
      end
   end
end

"""
    dump_data(io::IO, df::DataFrame)

Dump the raw data as text to `io`.
"""
function dump_data(io::IO, df::DataFrame)
   n = size(df,1)   
   for i=1:n
      outstr = @sprintf("%d %d %d %7.4f %d %d %d %d %d %d %10.4g %10.4g %10.4g",df.aid[i],df.sid[i],df.did[i],df.f[i],df.male[i],df.proven[i],df.dau[i],df.age[i],df.gen[i],df.nrecdau[i],df.ebv[i],df.ebv1st[i],df.ebvend[i])
      for k=1:length(df.bv[i])
         outstr = outstr * @sprintf(" %10.4g", ifelse(ismissing(df.bv[i][k]),0.0,df.bv[i][k]))
      end
      for k=1:length(df.pe[i])
         outstr = outstr * @sprintf(" %10.4g", ifelse(ismissing(df.pe[i][k]),0.0,df.pe[i][k]))
      end
      for k=1:length(df.te[i])
         outstr = outstr * @sprintf(" %10.4g", ifelse(ismissing(df.te[i][k]),0.0,df.te[i][k]))
      end
      for k=1:length(df.y[i])
         outstr = outstr * @sprintf(" %10.4g", ifelse(ismissing(df.y[i][k]),0.0,df.y[i][k]))
      end
      print(io,outstr,"\n")
   end
end

"""
    vg,ve = read_vc_1st(logfile; airemlf90=false)

Read genetic variance `vg` and residual variance `ve` from a log file generated by REMLF90.
Use `airemlf90=true` for a log file generated by AIREMLF90.
"""
function read_vc_1st(logfile; airemlf90=false)
   vg = 0
   ve = 0
   open(logfile,"r") do io
      readline(io)
      readline(io)
      if airemlf90
         readline(io)
         readline(io)
      end
      # vg
      readline(io)
      line = strip(readline(io))
      col = split(line)
      vg = parse(Float64,col[1])
      # ve
      readline(io)
      line = strip(readline(io))
      col = split(line)
      ve = parse(Float64,col[1])
   end
   return vg,ve
end

"""
    vg,vp,ve = read_vc_rep(logfile; airemlf90=false)

Read genetic variance `vg`, permanent environmental variance `vp`, and residual variance `ve` from a log file generated by REMLF90.
Use `airemlf90=true` for a log file generated by AIREMLF90.
"""
function read_vc_rep(logfile; airemlf90=false)
   vg = 0
   vp = 0
   ve = 0
   open(logfile,"r") do io
      readline(io)
      readline(io)
      if airemlf90
         readline(io)
         readline(io)
      end
      # vg
      readline(io)
      line = strip(readline(io))
      col = split(line)
      vg = parse(Float64,col[1])
      # vp
      readline(io)
      line = strip(readline(io))
      col = split(line)
      vp = parse(Float64,col[1])
      # ve
      readline(io)
      line = strip(readline(io))
      col = split(line)
      ve = parse(Float64,col[1])
   end
   return vg,vp,ve
end

"""
    G, E = read_vc_mt(logfile; airemlf90=false)

Read genetic (co)variance matrix `G` and residual (co)variance matrix `E` from a log file generated by REMLF90.
Use `airemlf90=true` for a log file generated by AIREMLF90.
"""
function read_vc_mt(logfile; airemlf90=false)
   ntr = 0
   open(logfile,"r") do io
      readline(io); readline(io)
      if airemlf90
         readline(io); readline(io)
      end
      # check the order of G
      readline(io)
      while(!eof(io))
         line = strip(readline(io))
         col = split(line)
         if isnothing(tryparse(Float64,col[1]))
            break
         end
         ntr = ntr + 1
      end
   end
   G = zeros(ntr,ntr)
   E = zeros(ntr,ntr)
   open(logfile,"r") do io
      readline(io)
      readline(io)
      if airemlf90
         readline(io)
         readline(io)
      end
      # G
      readline(io)
      load_matrix!(io,G)
      # skip some lines
      for i=1:(ntr*2 + 4)
         readline(io)
      end
      # E
      readline(io)
      load_matrix!(io,E)
   end
   return G,E
end

function load_matrix!(io,M)
   nr = size(M,1)
   nc = size(M,2)
   for i=1:nr
      line = strip(readline(io))
      col = split(line)
      for j=1:nc
         M[i,j] = parse(Float64,col[j])
      end
   end
end

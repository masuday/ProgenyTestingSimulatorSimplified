function write_data_rep(io::IO, df::DataFrame)
   nrec = size(df,1)
   nlac = length(df.y[1])
   for i=1:nrec
      for k=1:nlac
         if !ismissing(df.y[i][k])
            print(io, @sprintf("%d %d %.2f\n",i,1,df.y[i][k]))
         end
      end
   end
end

function write_data_mt(io::IO, df::DataFrame)
   nrec = size(df,1)
   nlac = length(df.y[1])
   for i=1:nrec
      nmiss = sum(ismissing.(df.y[i]))
      if nmiss < nlac
         print(io, @sprintf("%d %d",i,1))
         for k=1:nlac
            y = ifelse(ismissing(df.y[i][k]), 0.0, df.y[i][k])
            print(io, @sprintf(" %.2f",y))
         end
         print(io, @sprintf("\n"))
      end
   end
end

function write_pedigree(io::IO, df::DataFrame)
   nped = size(df,1)
   for i=1:nped
      aid = df.aid[i]
      sid = df.sid[i]
      did = df.did[i]
      inbcode = get_inbupg_code(sid,did,df.f)
      print(io, @sprintf("%d %d %d %4d\n",aid,sid,did,inbcode))
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

function write_parfile_rep(io::IO, df::DataFrame, datafile::String, pedfile::String, vg::Float64, vp::Float64, ve::Float64; option=Vector{String}[])
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
  2 1 cross
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
OPTION use_yams
OPTION EM-REML 5
")
   if length(option)>0
      for str in option
         print(io,str,"\n")
      end
   end
end

function write_parfile_mt(io::IO, df::DataFrame, datafile::String, pedfile::String, G::Matrix{Float64}, R::Matrix{Float64}; option=Vector{String}[])
   ntr = size(G,1)
   effect_line1 = repeat("1 ",ntr) * string(size(df,1)) * " cross"
   effect_line2 = repeat("2 ",ntr) * "1" * " cross"
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
   write_matrix(io,R)
   print(io,"RANDOM_GROUP
  1
RANDOM_TYPE
  add_an_upginb
FILE
  $(pedfile)
(CO)VARIANCES
")
   write_matrix(io,G)
   print(io,"OPTION use_yams
OPTION EM-REML 10
")

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
         print(io,@sprintf(" %.5f",M[i,j]))
      end
      print(io,@sprintf("\n"))
   end
end

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
         if lev>length(ebv)
            throw(DimensionMismatch("short in elements of ebv"))
         else
            ebv[trt,lev] = sol
         end
      end
   end
end


#
# empty dataframe
#
"""
   df = initial_data(gp::GeneticParameter)

Generate an initial (empty) daraframe `df` out of the genetic parameter `gp`.
"""
function initial_data(gp::GeneticParameter)
   ntr = gp.ntr
   df = DataFrame(
      aid    = Int64[],   # animal ID
      sid    = Int64[],   # sire ID
      did    = Int64[],   # dam ID
      f      = Float64[], # inbreeding coefficient
      herd   = Int64[],   # herd
      alive  = Bool[],    # alive or not
      male   = Bool[],    # male or not
      proven = Bool[],    # proven male or not
      dau    = Bool[],    # test daughter or not
      age    = Int64[],   # age
      gen    = Int64[],   # generation born
      nrecdau= Int64[],   # number of recorded daughters
      ebv    = Union{Missing,Float64}[],   # latest EBV
      ebv1st = Union{Missing,Float64}[],   # first-crop EBV
      ebv2nd = Union{Missing,Float64}[],   # second-crop EBV
      bv     = SVector{ntr,Float64}[],     # true BV
      pe     = SVector{ntr,Float64}[],     # PE
      te     = SVector{ntr,Float64}[],     # temporary environmental effect
      y      = Vector{Union{Float64,Missing}}[]  # phenotype
   )
   return df
end

#
# intial populations
#
"""
    generate_founders!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter; debug=false)

Generate founder males and females given genetic parameters `gp` and simulation parameters `sp`, and update the dataframe `df`.
The dataframe `df` must be initialized with `initial_data()`,
The debug option `debug` prints details of simulated data.
"""
function generate_founders!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter; debug=false)
   if debug; println("Generating founders"); end
   # bulls
   for age in reverse(1:sp.maxagem)
      if debug; println("  male age $(age): $(sp.nm[age]) bulls; ID $(size(df,1)+1) to $(size(df,1)+sp.nm[age])"); end
      for n=1:sp.nm[age]
         aid = size(df,1) + 1
         herd = 0
         bv,pe,te,y = generate_base_variables(gp)
         proven = ifelse(age>=agem_proven, true, false)
         # aid,sid,did,f,herd,alive,male,proven,dau,age,gen,nrecdau,ebv,ebv1st,ebv2nd,bv,pe,te,y
         push!(df, [aid,0,0,0.0,herd,true,true,proven,false,age,-age,0,missing,missing,missing,bv,pe,te,y])
      end
   end
   # cows
   for age in reverse(1:sp.maxagef)
      if debug; println("  female age $(age): $(sp.nf[age]) cows; ID $(size(df,1)+1) to $(size(df,1)+sp.nf[age])"); end
      for n=1:sp.nf[age]
         aid = size(df,1) + 1
         herd = 1
         bv,pe,te,y = generate_base_variables(gp)
         # aid,sid,did,f,herd,alive,male,proven,dau,age,gen,nrecdau,ebv,ebv1st,ebv2nd,bv,pe,te,y
         push!(df, [aid,0,0,0.0,herd,true,false,false,false,age,-age,0,missing,missing,missing,bv,pe,te,y])
      end
   end
   return nothing
end

#
# cow herd
#
"""
    assign_cow_herd!(df::DataFrame; method::Symbol, nherd=1)

This function assigns the herd code to each cow. The assignment method is defined by `method`, which may need extra arguments like `nherd`.
The following symbols can be allowed for the method.

- `:equal` (default): The cows are assigned to herds so that each hard has the (nearly) equal number of cows. It needs the extra argument `nherd`, the number of herds (default = 1).
- `:random`: The cows are randomly split into the herds. It also needs `nherd`.
"""
function assign_cow_herd!(df::DataFrame; method=:equal, nherd=1, debug=false)
   if method==:equal
      if debug; println("Assigning cows to $(nherd) herds (equally)"); end
      if nherd<1
         throw(ArgumentError("negative `nherd`"))
      end
      n = size(df,1)
      m = 1
      for i in 1:n
         if !df.male[i]
            df.herd[i] = m
            m = (m % nherd) + 1
         end
      end
   elseif method==:random
      if debug; println("Assigning cows to $(nherd) herds (randomly)"); end
      if nherd<1
         throw(ArgumentError("negative `nherd`"))
      end
      ncows = count(.!df.male)
      df.herd[.!df.male] .= rand(1:nherd,ncows)
   else
      throw(ArgumentError("unknown method $(method)"))
   end
end

function generate_base_variables(gp::GeneticParameter)
   bv = SVector{gp.ntr}(generate_random_vectors(gp.Lg))
   pe = SVector{gp.ntr}(generate_random_vectors(gp.Lp))
   te = SVector{gp.ntr}(generate_random_vectors(gp.Le))
   y = Vector{Union{Missing,Float64}}(fill(missing,gp.ntr))
   return bv,pe,te,y
end

function generate_bv(gp::GeneticParameter, sbv, dbv, fs, fd)
   ms = get_ms_std(fs,fd)
   bv = 0.5*sbv + 0.5*dbv + ms*( SVector{gp.ntr}(generate_random_vectors(gp.Lg)) )
   return bv
end

function get_ms_std(fs, fd)
   ms_std = sqrt( 0.5*(1 - 0.5*(fs+fd)) )
   return ms_std
end

function generate_random_vectors(L::Matrix{T}) where T
   n = size(L,1)
   return L*randn(n)
end

# animal index
function id_of_young_bulls(df::DataFrame)
   return df[df.alive .&& df.male .&& df.age.==agem_young, :aid]
end

function id_of_preselected_bulls(df::DataFrame)
   return df[df.alive .&& df.male .&& df.age.>=agem_young .&& df.age.<agem_proven, :aid]
end

function id_of_candidate_bulls(df::DataFrame)
   return df[df.alive .&& df.male .&& df.age.==agem_selection, :aid]
end

function id_of_proven_bulls(df::DataFrame)
   return df[df.alive .&& df.male .&& df.age.>=agem_proven .&& df.proven, :aid]
end

function id_of_cows(df::DataFrame)
   return df[df.alive .&& .!df.male .&& df.age.>=agef_first_lact, :aid]
end

function id_of_heifers_and_cows(df::DataFrame)
   return df[df.alive .&& .!df.male .&& df.age.>=agef_heifer, :aid]
end

function id_of_male_calves(df::DataFrame)
   return df[df.alive .&& df.male .&& df.age.==0, :aid]
end

function id_of_female_calves(df::DataFrame; which="all")
   if which=="test"
      return df[df.alive .&& .!df.male .&& df.age.==0 .&& df.dau, :aid]
   elseif which=="regular"
      return df[df.alive .&& .!df.male .&& df.age.==0 .&& .!df.dau, :aid]
   else
      return df[df.alive .&& .!df.male .&& df.age.==0, :aid]
   end
end

function cow_age_to_lactation(age)
   return age - agef_first_lact + 1
end

function select_top_id(key, id, n)
   return id[sortperm(key,rev=true)[1:n]]
end

"""
    list_of_top_young_bulls(df,n)

Returns a list of IDs for the top `n` bulls in each generation (year) given data `df`.
"""
function list_of_top_young_bulls(df,n)
   gyblist = Int[]
   for age in 1:3
      id   = df[df.alive .&& df.male .&& df.age.==age, :aid]
      gebv = df[df.alive .&& df.male .&& df.age.==age, :ebv1st]
      gyblist = [gyblist; id[sortperm(gebv,rev=true)[1:n]] ]
   end
   return unique(sort(gyblist))
end

"""
    function get_bull_index(df,bulls; list)

Returns bull indices from `df` conditioned by `bulls`, which should be a vector of the following symbols.

- `:young` = young bulls (age $(agem_young))
- `:preselected` = preselected bulls (age $(agem_young) - $(agem_selection))
- `:proven` = proven bulls (age >=$(agem_proven))
- `:custom` = bulls in `list`
"""
function get_bull_index(df,bulls; list=Int[])
   if !all( map(x->in(x,[:young,:preselected,:proven,:custom]), bulls) )
      throw( ArgumentError("Invalid argument") )
   end
   idx = Int[]
   if :young in bulls
      idx = [idx; id_of_young_bulls(df)]
   end
   if :preselected in bulls
      idx = [idx; id_of_preselected_bulls(df)]
   end
   if :proven in bulls
      idx = [idx; id_of_proven_bulls(df)]
   end
   if :custom in bulls
      idx = [idx; list]
   end
   return unique(sort(idx))
end

"""
    test_mating!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter, gen::Int64; herds=Integer[], bulls=[:young], debug=false)

Mate young bulls with randomly-selected cows to produce test daughters, and update the dataframe `df`.
The dataframe `df` should be prepared by `generate_founders!`.

Given `herds` as an integer array to specify a list of herds for females to be bred.

For `bulls`, see the docstrings for `ProgenyTestingSimulatorSimplified.get_bull_index`.
"""
function test_mating!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter, gen::Int64;
                      herds=Integer[], bulls=[:young], list=Int[], debug=false)
   if debug; println("Test mating in generation $(gen)"); end
   if debug; println("Bull category: $(bulls)"); end
   if debug && length(list)>0; println("Custom bull list: $(length(list))"); end
   # particular herds
   if length(herds)==0; herd_list=1:maximum(df.herd); else; herd_list=herds; end
   idx_testherd = map(x->issubset(x,herd_list), df.herd)
   # index for young bulls and cows
   # random mating
   #idx_ybulls = id_of_young_bulls(df)
   idx_ybulls = get_bull_index(df,bulls, list=list)
   idx_cows   = df[df.alive .&& .!df.male .&& idx_testherd .&& df.age.>=agef_first_lact, :aid]
   shuffle!(idx_cows)
   j = 0
   for sid in idx_ybulls
      for k=1:sp.ndau
         j = j + 1
         aid = size(df,1) + 1
         did = idx_cows[j]
         herd = df.herd[did]
         bv = generate_bv(gp, df.bv[sid], df.bv[did], df.f[sid], df.f[did])
         pe = SVector{gp.ntr}(generate_random_vectors(gp.Lp))
         te = SVector{gp.ntr}(generate_random_vectors(gp.Le))
         y = Vector{Union{Missing,Float64}}(fill(missing,gp.ntr))
         # aid,sid,did,f,herd,alive,male,proven,dau,age,gen,nrecdau,ebv,ebv1st,ebv2nd,bv,pe,te,y
         push!(df, [aid,sid,did,0.0,herd,true,false,false,true,0,gen,0,missing,missing,missing,bv,pe,te,y])
      end
   end
   if debug; println("  test $(j) daughters from $(length(idx_ybulls)) young and $(j) cows; ID $(size(df,1)-j+1) to $(size(df,1))"); end
   return nothing
end

"""
    regular_mating!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter, gen::Int64; sexratio=0.5, herds=Integer[], bulls=[:proven], debug=false)

Produce calves (male & female) by proven bulls with randomly-selected cows, and update the dataframe `df`.
The dataframe `df` should be prepared by `generate_founders!`.

Given `herds` as an integer array to specify a list of herds for females to be bred.

For `bulls`, see the docstrings for `ProgenyTestingSimulatorSimplified.get_bull_index`.
"""
function regular_mating!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter, gen::Int64;
                         sexratio=0.5, herds=Integer[], bulls=[:proven], list=Int[], debug=false)
   if sexratio<0.0 || sexratio>1
      throw(ArgumentError("sex ratio should be in range of 0<=x<=1"))
   end
   if debug; println("Regular mating in generation $(gen)"); end
   if debug; println("Bull category: $(bulls)"); end
   if debug && length(list)>0; println("Custom bull list: $(length(list))"); end
   # index for selected bulls and cows
   nm = 0
   nf = 0
   # random mating
   #idx_sbulls  = id_of_proven_bulls(df)
   idx_sbulls = get_bull_index(df,bulls, list=list)
   idx_allcows = df[df.alive .&& .!df.male .&& df.age.>=agef_heifer, :aid]
   idx_tested  = df[df.alive .&&   df.dau  .&& df.age.==0, :did]
   idx_cows    = setdiff(idx_allcows,idx_tested)
   if length(herds)==0; herd_list=1:maximum(df.herd); else; herd_list=herds; end
   for did in idx_cows
      dherd = df.herd[did]
      if dherd in herd_list
         aid = size(df,1) + 1
         sid = sample(idx_sbulls)
         male = rand()<sexratio ? true : false
         if male; nm=nm+1; else; nf=nf+1; end
         if male; herd=0; else; herd=df.herd[did]; end
         bv = generate_bv(gp, df.bv[sid], df.bv[did], df.f[sid], df.f[did])
         pe = SVector{gp.ntr}(generate_random_vectors(gp.Lp))
         te = SVector{gp.ntr}(generate_random_vectors(gp.Le))
         y = Vector{Union{Missing,Float64}}(fill(missing,gp.ntr))
         # aid,sid,did,f,herd,alive,male,proven,dau,age,gen,nrecdau,ebv,ebv1st,ebv2nd,bv,pe,te,y
         push!(df, [aid,sid,did,0.0,herd,true,male,false,false,0,gen,0,missing,missing,missing,bv,pe,te,y])
      end
   end
   if debug
      println("  total $(nm+nf) calves: $(nm) male and $(nf) female from $(length(idx_sbulls)) bulls and $(length(idx_cows)) cows; ID $(size(df,1)-length(idx_cows)+1) to $(size(df,1))")
   end
   return nothing
end

"""
    et_mating!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter, gen::Int64; sexratio=0.5, herds=Integer[], bulls=[:proven], dist=Returns(1), debug=false)

Produce multiple calves by mating proven bulls with randomly-selected cows using ET, and update the dataframe `df`.
The dataframe `df` should be prepared by `generate_founders!`.

Given `herds` as an integer array to specify a list of herds for females to be bred.

The function `dist` returns the number of embryos per cow (default = `Returns(1)`, always returns `1`).

For `bulls`, see the docstrings for `ProgenyTestingSimulatorSimplified.get_bull_index`.
"""
function et_mating!(df::DataFrame, gp::GeneticParameter, sp::SimulationParameter, gen::Int64; 
                    sexratio=0.5, herds=Integer[], bulls=[:proven], list=Int[], dist=Returns(1), debug=false)
   if sexratio<0.0 || sexratio>1
      throw(ArgumentError("sex ratio should be in range of 0<=x<=1"))
   end
   if debug; println("Regular ET mating in generation $(gen)"); end
   if debug; println("Bull category: $(bulls)"); end
   if debug && length(list)>0; println("Custom bull list: $(length(list))"); end
   # index for selected bulls and cows
   nm = 0
   nf = 0
   # random mating
   #idx_sbulls  = id_of_proven_bulls(df)
   idx_sbulls = get_bull_index(df,bulls, list=list)
   idx_allcows = df[df.alive .&& .!df.male .&& df.age.>=agef_heifer, :aid]
   idx_tested  = df[df.alive .&&   df.dau  .&& df.age.==0, :did]
   idx_cows    = setdiff(idx_allcows,idx_tested)
   if length(herds)==0; herd_list=1:maximum(df.herd); else; herd_list=herds; end
   for did in idx_cows
      dherd = df.herd[did]
      if dherd in herd_list
         aid = size(df,1) + 1
         sid = sample(idx_sbulls)
         nembryo = dist()
         for emb in 1:nembryo
            male = rand()<sexratio ? true : false
            if male; nm=nm+1; else; nf=nf+1; end
            if male; herd=0; else; herd=df.herd[did]; end
            bv = generate_bv(gp, df.bv[sid], df.bv[did], df.f[sid], df.f[did])
            pe = SVector{gp.ntr}(generate_random_vectors(gp.Lp))
            te = SVector{gp.ntr}(generate_random_vectors(gp.Le))
            y = Vector{Union{Missing,Float64}}(fill(missing,gp.ntr))
            # aid,sid,did,f,herd,alive,male,proven,dau,age,gen,nrecdau,ebv,ebv1st,ebv2nd,bv,pe,te,y
            push!(df, [aid,sid,did,0.0,herd,true,male,false,false,0,gen,0,missing,missing,missing,bv,pe,te,y])
         end
      end
   end
   if debug
      println("  total $(nm+nf) calves: $(nm) male and $(nf) female from $(length(idx_sbulls)) bulls and $(length(idx_cows)) cows; ID $(size(df,1)-length(idx_cows)+1) to $(size(df,1))")
   end
   return nothing
end

function get_pedigree_list(sires, dams; ml::Bool=false)
   if length(sires) != length(dams)
      throw(ArgumentError("unequal size of sire and dam vectors"))
   end
   n = length(sires)
   pedlist = zeros(Int,2,n)
   for i=1:n
      s = sires[i]
      d = dams[i]
      if ml==true
         pedlist[1,i] = max(s,d)
         pedlist[2,i] = min(s,d)
      else
         pedlist[1,i] = s
         pedlist[2,i] = d
      end
   end
   return pedlist
end

# kernel of Meuwissen and Luo (1992)
# specialized for new animals in the same generation
# the argument "ped" not being rewritten
function kernel_meuwissen_and_luo2!(ped::Matrix{Ti}, inb::Vector{Tv}; first=1) where {Tv<:Real, Ti<:Integer}
   n::Ti = size(ped,2)
   if n<length(inb)
      throw(DimensionMismatch("inb shorter than pedigree size"))
   end
   f = OffsetArray{Tv}(undef,0:n)
   f[1:n] .= inb[1:n]
   point = zeros(Ti,n,Threads.nthreads())
   T = zeros(Tv,n,Threads.nthreads())
   B = zeros(Tv,n)

   # start
   f[0] = -1
   for i=1:first-1
      s0 = ped[1,i]
      d0 = ped[2,i]
      B[i] = 0.5 - 0.25*(f[s0]+f[d0])
   end
   Base.Threads.@threads for i=first:n
      local s0 = ped[1,i]
      local d0 = ped[2,i]
      B[i] = 0.5 - 0.25*(f[s0]+f[d0])

      if s0==0 || d0==0
         f[i] = 0.0
      else
         local fi = -1.0
         T[i,Threads.threadid()] = 1.0
         local j = i
         while(j!=0)
            local k = j
            local r = 0.5 * T[k,Threads.threadid()]
            local ks = ped[1,k]
            local kd = ped[2,k]
            if ks != 0
               while(point[k,Threads.threadid()]>ks)
                  k = point[k,Threads.threadid()]
               end
               T[ks,Threads.threadid()] = T[ks,Threads.threadid()] + r
               if ks != point[k,Threads.threadid()]
                  point[ks,Threads.threadid()] = point[k,Threads.threadid()]
                  point[k,Threads.threadid()] = ks
               end
               if kd != 0
                  while(point[k,Threads.threadid()]>kd)
                     k = point[k,Threads.threadid()]
                  end
                  T[kd,Threads.threadid()] = T[kd,Threads.threadid()] + r
                  if kd != point[k,Threads.threadid()]
                     point[kd,Threads.threadid()] = point[k,Threads.threadid()]
                     point[k,Threads.threadid()] = kd
                  end
               end
            end
            fi = fi + T[j,Threads.threadid()]*T[j,Threads.threadid()]*B[j]
            T[j,Threads.threadid()] = 0.0
            k = j
            j = point[j,Threads.threadid()]
            point[k,Threads.threadid()] = 0
         end
         f[i] = fi
      end
   end
   inb[first:n] .= f[first:n]
   return nothing
end

# the latest generation only
"""
    update_inbreeding!(df)
    update_inbreeding!(df, age)

Update inbreeding for all animals in `df`. With age `age`, the computation will be limited to animals with the specified age.
"""
function update_inbreeding!(df::DataFrame, age::Int64=-1; debug=false)
   if debug; println("Inbreeding computation"); end
   pedlist = get_pedigree_list(df.sid, df.did, ml=true)
   if age<0
      firstidx = 1
   else
      firstidx = first(df[df.age.==age,:aid])
   end
   kernel_meuwissen_and_luo2!(pedlist, df.f, first=firstidx)
   if debug; println("  ID from $(firstidx) to $(size(df,1))"); end
   return nothing
end

# phenotyping
"""
    assign_phenotype!(df::DataFrame, gp::GeneticParameter; debug=false)

Assigning phenotype to cows by the model: *y = mu + bv + pe + te*.
"""
function assign_phenotype!(df::DataFrame, gp::GeneticParameter; debug=false)
   if debug; println("Assigning phenotype to cows"); end
   idx_cows = df[df.alive .&& .!df.male .&& df.age.>=agef_first_lact, :aid]
   n = 0
   for cow in idx_cows
      lact = df.age[cow] - agef_first_lact + 1
      # strict check: off
      #if !ismissing(df.y[cow][lact])
      #   throw(ErrorException("phenotype must have been missing."))
      #end
      n = n + 1
      df.y[cow][lact] = gp.mu[lact] + df.bv[cow][lact] + df.pe[cow][lact] + df.te[cow][lact]
      
      # mark at the first recording time only
      nownrec = sum(.!ismissing.(df.y[cow]))
      if nownrec==1
         sid = df.sid[cow]
         if sid>0
            df.nrecdau[sid] = df.nrecdau[sid] + 1
         end
         did = df.did[cow]
         if did>0
            df.nrecdau[did] = df.nrecdau[did] + 1
         end
      end
   end
   if debug; println("  generate phenotype for $(n) cows"); end
end


function show_current_population(df, sp)
   println("Current population structure")
   println(@sprintf("%4s %10s %10s %10s   %10s %10s %10s","","male","","","female","",""))
   println(@sprintf("%4s %10s %10s %10s   %10s %10s %10s","age","all","(proven)","[expected]","all","(testdau)","[expected]"))
   nm_all = length(df[df.alive .&& df.male .&& df.age.==0, :aid])
   nf_all = length(df[df.alive .&& .!df.male .&& df.age.==0, :aid])
   nf_dau = length(df[df.alive .&& df.dau .&& df.age.==0, :aid])
   ntestdau = expected_number_of_test_daughters(sp)
   ncalves = expected_number_of_calves(sp)
   nf_exp = Int(floor( 0.5*(ncalves-ntestdau) + ntestdau ))
   nm_exp = expected_number_of_calves(sp) - nf_exp
   println(@sprintf("%4d %10d (%8d) [%8d]   %10d (%8d) [%8d]",0,nm_all,NaN,nm_exp, nf_all,nf_dau,nf_exp))
   for age=1:max(sp.maxagem,sp.maxagef)
      nm_all    = length(df[df.alive .&& df.male .&& df.age.==age, :aid])
      nm_proven = length(df[df.alive .&& df.proven .&& df.age.==age, :aid])
      nf_all    = length(df[df.alive .&& .!df.male .&& df.age.==age, :aid])
      nf_dau    = length(df[df.alive .&& df.dau .&& df.age.==age, :aid])
      if age<=sp.maxagem; nm_exp=sp.nm[age]; else nm_exp=0; end
      if age<=sp.maxagef; nf_exp=sp.nf[age]; else nf_exp=0; end
      println(@sprintf("%4d %10d (%8d) [%8d]   %10d (%8d) [%8d]",age,nm_all,nm_proven,nm_exp, nf_all,nf_dau,nf_exp))
   end
end

function expected_number_of_calves(sp)
   nyb = sp.nm[agem_young]
   ntestdau = sp.ndau * nyb
   # possible cows for test daughters
   nallcows = sum(sp.nf[agef_first_lact:end])
   #nothercows = nallcows - ntestdau
   return nallcows + sp.nf[agef_heifer]
end

function expected_number_of_test_daughters(sp)
   nyb = sp.nm[agem_young]
   return sp.ndau * nyb
end

function mean_vector(v)
   avg = zeros(length(v[1]))
   for i=1:length(v)
      avg .= avg .+ v[i]
   end
   avg .= avg / length(v)
   return avg
end

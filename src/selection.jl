"""
    selected_id = select_animals(id,crit,nsel::Int64; desc=true)

Return top `nsel` IDs out of animal ID `id` based on criterion `crit`.
By default, the criterion will be sorted by descending order (high to low).
"""
function select_animals(id,crit,nsel; desc=true)
   if length(id)!=length(crit)
      throw(DimensionMismatch("inconsistent length of animal id and criterion"))
   end
   selected_id = id[ sortperm(crit,rev=desc)[1:min(nsel,length(id))] ]
   return selected_id
end

function candidate_bull_selection!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)
   id = id_of_candidate_bulls(df)
   nall = length(id)
   if debug; println("Selection of candidates: $(nall) bulls"); end
   if random
      if debug; println("  random selection"); end
      ebv = rand(nall)
   else
      if debug; println("  by EBV"); end
      ebv = df.ebv[id]
   end
   perm = sortperm(ebv,rev=true)
   # keep the first nsel=sp.nm[agem_proven] bulls;
   nsel = sp.nm[agem_proven]
   # cull the remaining bulls.
   for i=nsel+1:nall
      df.alive[ id[perm[i]] ] = false
   end
   if debug; println("  culled $(nall-nsel) bulls"); end
   nothing
end

function male_calf_selection!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)
   id = id_of_male_calves(df)
   nall = length(id)
   if debug; println("Selection of male calves: $(nall) calves in total"); end
   if random
      if debug; println("  random selection"); end
      ebv = rand(nall)
   else
      if debug; println("  by EBV"); end
      ebv = df.ebv[id]
   end
   perm = sortperm(ebv,rev=true)
   # keep the first nsel=sp.nm[agem_young] bulls;
   nsel = sp.nm[agem_young]
   # cull the remaining bulls.
   n = 0
   for i=nsel+1:nall
      n = n + 1
      df.alive[ id[perm[i]] ] = false
   end
   if debug; println("  culled $(n) male calves - $(nall-n) calves survive"); end
   nothing
end

function female_calf_selection!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)
   # test daughters must not be culled; only regular daughter calves
   id = id_of_female_calves(df, which="regular")
   nall = length(id)
   if debug; println("Selection of female calves: $(nall) calves in total"); end
   if random
      if debug; println("  random selection"); end
      ebv = rand(nall)
   else
      if debug; println("  by EBV"); end
      ebv = df.ebv[id]
   end
   perm = sortperm(ebv,rev=true)
   # keep the first nsel=sp.nf[agef_heifer] animals;
   nsel = sp.nf[agef_heifer]
   # cull the remaining bulls.
   n = 0
   for i=nsel+1:nall
      n = n + 1
      df.alive[ id[perm[i]] ] = false
   end
   if debug; println("  culled $(n) female calves - $(nall-n) calves survive"); end
end

function cull_old_bulls!(df::DataFrame, sp::SimulationParameter; debug=false)
   if debug; println("Culling the oldest bulls"); end
   id = df[df.alive .&& df.male .&& df.age.==sp.maxagem, :aid]
   df[id,:alive] .= false
   if debug; println("  age $(sp.maxagem): culled $(length(id)) bulls, all"); end
   nothing
end

function cull_some_heifers_and_cows!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)
   if debug; println("Culling some cows/heifers:"); end
   #for age=agef_first_lact:sp.maxagef
   for age=agef_heifer:sp.maxagef
      if age==agef_heifer
         # exclude test daughters
         status = "heifers"
         id = df[df.alive .&& .!df.male .&& df.age.==age .&& .!df.dau, :aid]
      else
         # for all cows
         status = "cows"
         id = df[df.alive .&& .!df.male .&& df.age.==age, :aid]
      end
      nall = length(id)
      if age==sp.maxagef
         # cull all the oldest cows
         df[id,:alive] .= false
         if debug; println("  age $(age): culled $(length(id)) cows, all"); end
      else
         # remove some cows
         nsel = min(sp.nf[age+1], nall)
         if random
            method = "by random"
            ebv = rand(nall)
         else
            method = "by EBV"
            ebv = df.ebv[id]
         end
         perm = sortperm(ebv,rev=true)
         # keep the first nsel cows;
         # cull the remaining cows.
         ncull = 0
         for i=nsel+1:nall
            ncull = ncull + 1
            df.alive[ id[perm[i]] ] = false
         end
         if debug; println("  age $(age): culled $(ncull) " * status * " " * method); end
      end
   end
   nothing
end

function increment_age!(df::DataFrame)
   id = df[df.alive, :aid]
   df.age[id] .= df.age[id] .+ 1
   nothing
end

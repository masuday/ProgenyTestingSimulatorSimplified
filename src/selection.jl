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

"""
    candidate_bull_selection!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)

Select candidate bulls (at age 4) and update `df`.
By default, `random=false`, the candidates will be selected by EBV; with `random=true`, random sampling will be applied.
"""
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
   for i=1:nsel
      df.proven[ id[perm[i]] ] = true
   end
   # cull the remaining bulls.
   for i=nsel+1:nall
      df.alive[ id[perm[i]] ] = false
   end
   if debug; println("  culled $(nall-nsel) bulls"); end
   nothing
end

"""
    male_calf_selection!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)

Select male calves (at age 0) and update `df`.
By default, `random=false`, the calves will be selected by EBV (i.e., PA); with `random=true`, random sampling will be applied.
"""
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

"""
    female_calf_selection!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)

Select female calves (at age 0) and update `df`.
By default, `random=false`, the calves will be selected by EBV (i.e., PA); with `random=true`, random sampling will be applied.
"""
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

"""
    cull_old_bulls!(df::DataFrame, sp::SimulationParameter; debug=false)

Cull the oldest, proven bulls (age=7) and update `df`.
"""
function cull_old_bulls!(df::DataFrame, sp::SimulationParameter; debug=false)
   if debug; println("Culling the oldest bulls"); end
   id = df[df.alive .&& df.male .&& df.age.==sp.maxagem, :aid]
   df[id,:alive] .= false
   if debug; println("  age $(sp.maxagem): culled $(length(id)) bulls, all"); end
   nothing
end

"""
    cull_some_heifers_and_cows!(df::DataFrame, sp::SimulationParameter; random=false, debug=false)

Cull heifers and cows (age=1 to 6) and update `df`.
The number of selected animals is defined in `sp`.
By default, `random=false`, the calves will be selected by EBV (i.e., PA); with `random=true`, random sampling will be applied.
"""
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

"""
    increment_age!(df::DataFrame)

Add 1 to age for all animals alive.
"""
function increment_age!(df::DataFrame)
   id = df[df.alive, :aid]
   df.age[id] .= df.age[id] .+ 1
   nothing
end

"""
    save_first_crop_ebv!(df::DataFrame)

Save the first-crop EBV for proven bulls and recorded cows.
"""
function save_first_crop_ebv!(df::DataFrame)
   n = size(df,1)
   for i=1:n
      if df.male[i]
         # bulls with recorded daughters
         if ismissing(df.ebv1st[i]) && df.nrecdau[i]>0
            df.ebv1st[i] = df.ebv[i]
         end
      else
         # recorded female
         if ismissing(df.ebv1st[i]) && sum(.!ismissing.(df.y[i]))>=1
            df.ebv1st[i] = df.ebv[i]
         end
      end
   end
end

"""
    save_second_crop_ebv!(df::DataFrame)

Save the second-crop EBV for proven bulls and recorded cows.
"""
function save_second_crop_ebv!(df::DataFrame, sp::SimulationParameter)
   n = size(df,1)
   for i=1:n
      if df.male[i]
         # bulls with recorded daughters (more than ndau after progeny testing)
         if ismissing(df.ebv2nd[i]) && df.nrecdau[i]>sp.ndau
            df.ebv2nd[i] = df.ebv[i]
         end
      else
         # recorded female (2 or more records)
         if ismissing(df.ebv2nd[i]) && sum(.!ismissing.(df.y[i]))>=2
            df.ebv2nd[i] = df.ebv[i]
         end
      end
   end
end

"""
    drop_culled_calves!(df::DataFrame)

Remove culled calves (age=0) and recode ID for selected calves.
This function must be called before `increment_age!`.
"""
function drop_culled_calves!(df::DataFrame; debug=false)
   if debug; println("Removing culled calves from the database"); end
   firstrow = findfirst(df.age.==0)
   nalive = 0
   if !isnothing(firstrow)
      # remove all culled calves from the dataframe
      culled_id = df[.!df.alive .&& df.age.==0,:aid]
      filter!(:aid=>x->!(x in culled_id), df)
      if debug; println("  $(length(culled_id)) calves dropped"); end
      # recode ID for the remaining calves
      for i=firstrow:size(df,1)
         df[i,:aid] = firstrow + nalive
         nalive = nalive + 1
      end
   end
   if debug; println("  $(nalive) calves survived"); end
   return nothing
end


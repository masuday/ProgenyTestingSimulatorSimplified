using ProgenyTestingSimulatorSimplified
using LinearAlgebra
using Test

function test_parameters()
   ntr = 5
   G = [
      10  8  6  4  2
       8 10  8  6  4
       6  8 10  8  6
       4  6  8 10  8
       2  4  6  8 10
   ] * 1.0
   P = copy(G)
   E = copy(G)
   gp = GeneticParameter(G,P,E)

   nm = [  4,  4,  4,  2,  2,  2]
   nf = [  6,  6,  4,  3,  2,  1]
   maxlact = ntr
   ndau = 2
   sp = SimulationParameter(nm,nf,maxlact,ndau)

   # aid,sid,did,f,alive,male,proven,dau,age,gen,ebv,ebv1st,bv,pe,te,y
   #df = initial_data(gp)
   #generate_founders!(df,gp,sp, debug=true)

   return gp,sp
end

function generate_blup_files()
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   test_mating!(df,gp,sp,0, debug=false)
   regular_mating!(df,gp,sp,0, debug=false)
   assign_phenotype!(df, gp, debug=false)
   
   open("ped.txt","w") do io
      write_pedigree(io,df)
   end
   open("rep.txt","w") do io
      write_data_rep(io,df)
   end
   open("mt.txt","w") do io
      write_data_mt(io,df)
   end
   open("rep.par","w") do io
      write_parfile_rep(io, df, "rep.txt", "ped.txt", 3.0, 3.0, 3.0)
   end
   open("mt.par","w") do io
      write_parfile_mt(io, df, "mt.txt", "ped.txt", gp.G, gp.P+gp.E)
   end
   run(`blupf90 rep.par`)
   run(`blupf90 mt.par`)
end

@testset "founder population" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   # male
   @test all( df[1:2,:age] .== 6 )
   @test all( df[3:4,:age] .== 5 )
   @test all( df[5:6,:age] .== 4 )
   @test all( df[7:10,:age] .== 3 )
   @test all( df[11:14,:age] .== 2 )
   @test all( df[15:18,:age] .== 1 )

   @test all( df[1:2,:gen] .== -6 )
   @test all( df[3:4,:gen] .== -5 )
   @test all( df[5:6,:gen] .== -4 )
   @test all( df[7:10,:gen] .== -3 )
   @test all( df[11:14,:gen] .== -2 )
   @test all( df[15:18,:gen] .== -1 )
   @test all( df[1:18,:male] )
   # female
   @test all( df[19:19,:age] .== 6 )
   @test all( df[20:21,:age] .== 5 )
   @test all( df[22:24,:age] .== 4 )
   @test all( df[25:28,:age] .== 3 )
   @test all( df[29:34,:age] .== 2 )
   @test all( df[35:40,:age] .== 1 )
   @test all( df[19:19,:gen] .== -6 )
   @test all( df[20:21,:gen] .== -5 )
   @test all( df[22:24,:gen] .== -4 )
   @test all( df[25:28,:gen] .== -3 )
   @test all( df[29:34,:gen] .== -2 )
   @test all( df[35:40,:gen] .== -1 )
   @test all( .!df[19:40,:male] )
   # variables
   for i=1:40
      @test length(df[i,:bv])==5
      @test length(df[i,:pe])==5
      @test length(df[i,:te])==5
      @test length(df[i,:y])==5
      @test all( ismissing.(df[i,:y]) )
   end
end

@testset "mating" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)

   # test mating: the number of test daughters: 4 bulls * 2 each = 8
   test_mating!(df,gp,sp,0, debug=false)
   @test size(df,1)==40+8
   @test all( df[41:48,:dau] .== true )
   @test all( df[41:48,:male] .== false )
   @test all( df[41:48,:age] .== 0 )
   @test all( df[41:48,:gen] .== 0 )
   @test all( df[41:48,:sid] .== [15,15,16,16,17,17,18,18])
   @test issubset(df[41:48,:did],19:40)

   # regular mating: the number of daughters: 6+(6+4+3+2+1)-8 = 14
   regular_mating!(df,gp,sp,0, debug=false)
   @test size(df,1)==48+6+(6+4+3+2+1)-8
   @test issubset(df[49:62,:sid],1:6)
   @test issubset(df[49:62,:did],setdiff(19:40,df[41:18,:did]))
   @test all( df[49:62,:dau] .== false )
   @test all( df[49:62,:age] .== 0 )
   @test all( df[49:62,:gen] .== 0 )
end

@testset "animal ID" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   test_mating!(df,gp,sp,0, debug=false)
   regular_mating!(df,gp,sp,0, debug=false)

   @test issubset(ProgenyTestingSimulatorSimplified.id_of_young_bulls(df), 15:18)
   @test issubset(ProgenyTestingSimulatorSimplified.id_of_candidate_bulls(df), 7:10)
   @test issubset(ProgenyTestingSimulatorSimplified.id_of_proven_bulls(df), 1:6)
   @test issubset(ProgenyTestingSimulatorSimplified.id_of_cows(df), 19:34)
   @test issubset(ProgenyTestingSimulatorSimplified.id_of_heifers_and_cows(df), 19:40)
   @test issubset(ProgenyTestingSimulatorSimplified.id_of_male_calves(df), df[df.alive.==true .&& df.age.==0 .&& df.male,:aid])
   @test issubset(ProgenyTestingSimulatorSimplified.id_of_female_calves(df), df[df.alive.==true .&& df.age.==0 .&& .!df.male,:aid])
end

@testset "phenotype" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   test_mating!(df,gp,sp,0, debug=false)
   regular_mating!(df,gp,sp,0, debug=false)
   assign_phenotype!(df, gp, debug=false)

   id_cows = ProgenyTestingSimulatorSimplified.id_of_cows(df)
   for aid=1:62
      if aid in id_cows
         lact = ProgenyTestingSimulatorSimplified.cow_age_to_lactation(df.age[aid])
         @test !ismissing(df.y[aid][lact])
      else
         @test sum(ismissing.(df.y[aid]))==5
      end
   end
end

@testset "first-crop daughters" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   for gen=0:4
      test_mating!(df,gp,sp,gen, debug=false)
      regular_mating!(df,gp,sp,gen, debug=false)
      update_inbreeding!(df,0, debug=false)
      assign_phenotype!(df,gp, debug=false)

      # random selection
      candidate_bull_selection!(df, sp, random=true, debug=false)
      male_calf_selection!(df, sp, random=true, debug=false)
      female_calf_selection!(df, sp, random=true, debug=false)

      # random culling
      cull_old_bulls!(df,sp, debug=false)
      cull_some_heifers_and_cows!(df, sp, random=true, debug=false)
      drop_culled_calves!(df, debug=false)
      increment_age!(df)
   end

   # confirm candidate bulls have test daughters with records in the 1st lactation
   # candidate bulls born in year 1: test daughters born in year 2; first calving in year 4
   bull_id = df[df.gen.==1 .&& df.male,:aid]
   dau_id = findall(x->x in bull_id, df.sid)
   @test all(df[dau_id,:dau])
   @test all(df[dau_id,:gen] .== 2)
   for i in dau_id
      @test sum(.!ismissing.(df[i,:y])) == 1
   end
end

@testset "variance components" begin
   # simple animal model
   (vg,ve) = read_vc_1st("remlf90.1st")
   @test vg ≈ 43.13
   @test ve ≈ 71.34
   (vg,ve) = read_vc_1st("airemlf90.1st", airemlf90=true)
   @test vg ≈ 43.081
   @test ve ≈ 71.378
   # repeatability model
   (vg,vp,ve) = read_vc_rep("remlf90.rep")
   @test vg ≈ 50.72
   @test vp ≈ 26.14
   @test ve ≈ 79.00
   (vg,vp,ve) = read_vc_rep("airemlf90.rep", airemlf90=true)
   @test vg ≈ 50.821
   @test vp ≈ 26.046
   @test ve ≈ 79.006
   # MT model
   refG = [
      43.894       41.387       39.208       37.167       35.293    
      41.387       44.129       45.256       44.994       43.390    
      39.208       45.256       48.514       49.464       48.167    
      37.167       44.994       49.464       51.205       50.272    
      35.293       43.390       48.167       50.272       49.758    
   ]
   refE = [
      70.700       25.810       23.637       21.915       14.582    
      25.810       106.97       35.393       35.011       35.773    
      23.637       35.393       123.97       47.822       50.766    
      21.915       35.011       47.822       136.46       58.051    
      14.582       35.773       50.766       58.051       136.47    
   ]
   (G,E) = read_vc_mt("remlf90.mt")
   @test refG ≈ G
   @test refE ≈ E
   (G,E) = read_vc_mt("airemlf90.mt", airemlf90=true)
   @test refG ≈ G
   @test refE ≈ E
end

@testset "solutions" begin
   # repeatability model
   ebv_ref  = [1.11, 2.22, 3.33, 4.44, 5.55]
   ebv_curr = zeros(5)
   load_solutions_rep!("test_solutions.rep", ebv_curr)
   @test isapprox(ebv_ref,ebv_curr, rtol=1e-4)

   @test_throws ErrorException load_solutions_rep!("test_solutions.mt", ebv_curr)

   # multiple trait model
   ebv_ref  = [1.11 2.11 3.11 4.11 5.11
               1.22 2.22 3.22 4.22 5.22]
   ebv_curr = zeros(2,5)
   load_solutions_mt!("test_solutions.mt", ebv_curr)
   @test isapprox(ebv_ref,ebv_curr, rtol=1e-4)

   # multiple trait model with limited traits
   ebv_ref  = [1.11 2.11 3.11 4.11 5.11]
   ebv_curr = zeros(1,5)
   load_solutions_mt!("test_solutions.st", ebv_curr)
   @test isapprox(ebv_ref,ebv_curr, rtol=1e-4)
end

@testset "running BLUP" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   test_mating!(df,gp,sp,0, debug=false)
   regular_mating!(df,gp,sp,0, debug=false)
   assign_phenotype!(df, gp, debug=false)

   # repeatability model
   #run(`blupf90 rep.par`)
   run(pipeline(`blupf90 rep.par`,devnull))
   n = size(df,1)
   ebv_ref = zeros(n)
   load_solutions_rep!("solutions.rep", ebv_ref)
   ebv_curr = zeros(n)
   load_solutions_rep!("solutions", ebv_curr)
   @test isapprox(ebv_ref,ebv_curr, rtol=1e-4)
   rm("solutions")

   # MT model
   #run(`blupf90 mt.par`)
   run(pipeline(`blupf90 mt.par`,devnull))
   n = size(df,1)
   ebv_ref = zeros(5,n)
   load_solutions_mt!("solutions.mt", ebv_ref)
   ebv_curr = zeros(5,n)
   load_solutions_mt!("solutions", ebv_curr)
   @test isapprox(ebv_ref,ebv_curr, rtol=1e-4)
   rm("solutions")
end

@testset "basic selection" begin
   id  = [ 1, 3, 5, 7, 9,11,13]
   ebv = [12,14,15,17,16,13,11]*1.0
   selected_id = ProgenyTestingSimulatorSimplified.select_animals(id,ebv,5)
   @test all( selected_id .== [7,9,5,3,11])
end

@testset "selection" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   test_mating!(df,gp,sp,0, debug=false)
   regular_mating!(df,gp,sp,0, debug=false)
   assign_phenotype!(df, gp, debug=false)
   n = size(df,1)
   df.ebv = rand(n)

   # candidate bulls: 
   id = ProgenyTestingSimulatorSimplified.id_of_candidate_bulls(df)
   p = sortperm(df.ebv[id],rev=true)
   selected = id[p[1:2]]
   nonselected = id[p[3:4]]
   candidate_bull_selection!(df, sp, random=false, debug=false)
   @test all( df.alive[selected] )
   @test all( .!(df.alive[nonselected]) )

   # male calves
   id = ProgenyTestingSimulatorSimplified.id_of_male_calves(df)
   p = sortperm(df.ebv[id],rev=true)
   selected = id[p[1:4]]
   nonselected = id[p[5:end]]
   male_calf_selection!(df, sp, random=false, debug=false)
   @test all( df.alive[selected] )
   @test all( .!(df.alive[nonselected]) )

   # female calves
   id_all = ProgenyTestingSimulatorSimplified.id_of_female_calves(df)
   id_dau = df[df.alive .&& .!df.male .&& df.dau .&& df.age.==0, :aid]
   id = setdiff(id_all,id_dau)
   p = sortperm(df.ebv[id],rev=true)
   ncalves = min(length(id),6)
   selected = id[p[1:ncalves]]
   nonselected = id[p[ncalves+1:end]]
   female_calf_selection!(df, sp, random=false, debug=false)
   @test all( df.alive[selected] )
   if length(nonselected)>0
      @test all( .!(df.alive[nonselected]) )
   end
end

@testset "culling" begin
   gp,sp = test_parameters()
   df = initial_data(gp)
   generate_founders!(df,gp,sp, debug=false)
   test_mating!(df,gp,sp,0, debug=false)
   regular_mating!(df,gp,sp,0, debug=false)
   assign_phenotype!(df, gp, debug=false)
   n = size(df,1)
   df.ebv = rand(n)

   # cull old bulls
   cull_old_bulls!(df,sp)
   @test( all( .!df[1:2,:alive]) )
   @test( all( df[3:end,:alive]) )

   # cull cows
   cull_some_heifers_and_cows!(df, sp, random=false, debug=false)
   # heifers
   id = df[df.alive .&& .!df.male .&& df.age.==1, :aid]
   @test sum(df[id,:alive])==sp.nf[2]
   # 1st lactation cows
   id = df[df.alive .&& .!df.male .&& df.age.==2, :aid]
   @test sum(df[id,:alive])==sp.nf[3]
   # 2nd lactation cows
   id = df[df.alive .&& .!df.male .&& df.age.==3, :aid]
   @test sum(df[id,:alive])==sp.nf[4]
   # 3rd lactation cows
   id = df[df.alive .&& .!df.male .&& df.age.==4, :aid]
   @test sum(df[id,:alive])==sp.nf[5]
   # 4th lactation cows
   id = df[df.alive .&& .!df.male .&& df.age.==5, :aid]
   @test sum(df[id,:alive])==sp.nf[6]
   # 5th lactation cows
   id = df[df.alive .&& .!df.male .&& df.age.==6, :aid]
   @test sum(df[id,:alive])==0
end

@testset "fix matrix by eigenvalues" begin
   G = [
      43.002  41.387  39.208  37.167  35.293
      41.387  44.129  45.256  44.994  43.39
      39.208  45.256  48.514  49.464  48.167
      37.167  44.994  49.464  51.205  50.272
      35.293  43.39   48.167  50.272  49.758
   ]
   fixedG = copy(G)
   fix_covariance_matrix!(fixedG)
   @test all(eigen(fixedG).values .> 0.0)
   @test isapprox(G,fixedG, rtol=0.001)
end

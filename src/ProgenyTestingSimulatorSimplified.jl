module ProgenyTestingSimulatorSimplified

using Base.Threads
using Random
using Printf
using LinearAlgebra
using StatsBase
using DataFrames
using OffsetArrays
using StaticArrays
using PedigreeBase
using SparseMatrixDicts

const agem_young       = 1
const agem_selection   = 4
const agem_proven      = 5
const agef_heifer      = 1
const agef_first_lact  = 2

export GeneticParameter, SimulationParameter
export initial_data, generate_founders!, test_mating!, regular_mating!, update_inbreeding!, assign_phenotype!
export write_data_1st, write_data_rep, write_data_mt, write_pedigree, 
   write_parfile_1st, write_parfile_rep, write_parfile_mt, 
   load_solutions_1st!, load_solutions_rep!, load_solutions_mt!, dump_data,
   read_vc_1st, read_vc_rep, read_vc_mt
export candidate_bull_selection!, male_calf_selection!, female_calf_selection!,
   cull_old_bulls!, cull_some_heifers_and_cows!, increment_age!, 
   save_first_crop_ebv!, save_second_crop_ebv!, drop_culled_calves!

#
# structure for genetic parameters
#
struct GeneticParameter{T<:Real}
   ntr::Int64
   G::Matrix{T}
   P::Matrix{T}
   E::Matrix{T}
   Lg::Matrix{T}
   Lp::Matrix{T}
   Le::Matrix{T}
   mu::Vector{T}
   function GeneticParameter(G::Matrix{T},P::Matrix{T},E::Matrix{T}) where T
      ntr = size(G,1)
      Lg = Matrix(cholesky(G).L)
      Lp = Matrix(cholesky(P).L)
      Le = Matrix(cholesky(E).L)
      mu = ones(ntr) .* 50.0
      new{T}(ntr,G,P,E,Lg,Lp,Le,mu)
   end
end

#
# simulation parameters
#
struct SimulationParameter
   nm::Vector{Int64} # number of male at each age
   nf::Vector{Int64} # number of female at each age
   maxagem::Int64    # max. age of male
   maxagef::Int64    # max. age of female
   maxlact::Int64    # number of lactations
   ndau::Int64       # number of tested daughters per young bull
   npb::Int64        # number of proven bulls
   function SimulationParameter(nm,nf,maxlact,ndau)
      maxagem = length(nm)
      maxagef = length(nf)
      if maxagem < agem_proven
         throw(DimensionMismatch("Length of nm must be >=$(agem_proven)."))
      end
      if maxagef < agef_first_lact
         throw(DimensionMismatch("Length of nf must be >=$(agef_first_lact)."))
      end
      npb = sum(nm[agem_proven:end])
      new(nm,nf,maxagem,maxagef,maxlact,ndau,npb)
   end
end

include("data.jl")
include("io.jl")
include("selection.jl")
include("tools.jl")

end

using JWAS,JWAS.Datasets,DataFrames,CSV,LinearAlgebra,StatsBase,Statistics,Random,Distributions,DelimitedFiles

QTL_number = [30,90,300];
QTL_index = 1  # QTL_index = 1, 2, 3
nqtl = QTL_number[QTL_index]

replicate = 1  # replicate = 1 to 30, 30 replicates for each QTL number

geno = CSV.read("controlled_genotype.csv", delim = ',', header=true, missingstrings=["NA"])
geno = geno[2:end]
nInds, nMarkers = size(geno)  # 413 individuals with 33701 markers

chr_dense = CSV.read("dense_pos.csv", delim = ',', header=true)  # dense windows were determined first
QTL_space = filter!(x->x∉chr_dense[:,2],Vector{Int64}(1:nMarkers))  # simulate QTL effects on nondense windows
Random.seed!(QTL_index * 100 + replicate)
QTL_pos = sample(QTL_space, nqtl, replace = false)
sort!(QTL_pos)

## QTL effects simulation from a standard normal distribution
Random.seed!(QTL_index * 100 + replicate)
trait_effect0 = randn(nqtl, 1)
G = 1
GU = cholesky(G).U
trait_effect = trait_effect0 * GU
QTL_effect = trait_effect

## Breeding value simulation
BV0 = Matrix(geno)[:,QTL_pos] * QTL_effect
std_matrix = Diagonal(vec(std(BV0,dims = 1)))   # make genetic variance = 1 for convenience
BV = BV0 * inv(std_matrix)

## Phenotype simulation
Random.seed!(QTL_index * 100 + replicate)
res = randn(nInds)
h2 = 0.5
sigma_g2 = 1   # genetic variance is 1 by standardizing BV
vare = sigma_g2 * (1 - h2)/h2
RU = cholesky(vare).U
residual = res * RU
pheno = BV + residual

Pheno_data = DataFrame(Phenotype = pheno[:,1], BV = BV[:,1], residual = residual[:,1])

QTL_data = DataFrame(QTL_pos = QTL_pos, QTL_effect = QTL_effect[:,1])

work_path = "/BayesC_90_simulation/"

mkpath(work_path*"nqtl_"*string(QTL_number[QTL_index])*"_replicate_"*string(replicate))

result_path = work_path*"nqtl_"*string(nqtl)*"_replicate_"*string(replicate)*"/";

phenotypes = DataFrame(NSFTVID = 1:nrow(Pheno_data), y = pheno[:,1]);

cd(result_path)

model_equation ="y = intercept"

R = 1.0
model = build_model(model_equation, R);

G3 = 1.0
add_genotypes(model, work_path * "controlled_genotype.csv",G3,separator=',');

Pi = 0.99
out_sem=runMCMC(model,phenotypes,methods="BayesC",Pi=Pi,estimatePi=true,chain_length=50000,
    output_samples_frequency=1000,burnin = 10000);

## GWAS using constant T (threshold = 1 / number of windows)
GWAS1 = GWAS(model,work_path*"controlled_map.csv",result_path*"MCMC_samples_marker_effects_y.txt",
    header=true,window_size="1 Mb",GWAS=true,threshold=1/378)
CSV.write("window_1MB.csv", GWAS1[1])

## GWAS using window specific T
function GWAS_mod(mme,map_file::AbstractString,marker_effects_file::AbstractString...;
              #window
              window_size = "1 Mb",sliding_window = false,
              #GWAS
              GWAS = true, threshold = 0.001,
              #genetic correlation
              genetic_correlation = false,
              #misc
              header = true, output_winVarProps = false)

    if split(window_size)[2] != "Mb"
        error("The format for window_size is \"1 Mb\".")
    end
    if map_file == false
        println("The map file is not provided. A fake map file is generated with 100 markers in each 1 Mb window.")
        nmarkers=length(readdlm(marker_effect_file,',',header=true)[2])
        mapfile = DataFrame(markerID=1:nmarkers,
                            chromosome=fill(1,nmarkers),
                            position=1:10_000:nmarkers*10_000)
        CSV.write("mapfile.temp",mapfile)
    end

    window_size_Mb = map(Int64,parse(Float64,split(window_size)[1])*1_000_000)
    mapfile = (header == true ? readdlm(map_file,',',header=true)[1] : readdlm(map_file,','))
    chr     = map(string,mapfile[:,2])
    pos     = map(Int64,mapfile[:,3])

    window_size_nSNPs   = Array{Int64,1}()  #save number of markers in ith window for all windows
    window_chr          = Array{String,1}() #1
    window_pos_start    = Array{Int64,1}()  #1_000_000
    window_pos_end      = Array{Int64,1}()  #2_000_000
    window_snp_start    = Array{Int64,1}()  #1_314_314
    window_snp_end      = Array{Int64,1}()  #1_999_003
    window_column_start = Array{Int64,1}()  #101
    window_column_end   = Array{Int64,1}()  #200

    index_start = 1
    for i in unique(chr)
      pos_on_chri     = pos[chr.== i] #assume chr and pos are sorted
      if sliding_window == false
          nwindow_on_chri = ceil(Int64,pos_on_chri[end]/window_size_Mb)
      else
          nwindow_on_chri = findfirst(x -> x >= pos_on_chri[end] - window_size_Mb, pos_on_chri)
      end

      for j in 1: nwindow_on_chri
        if sliding_window == false
            thisstart = window_size_Mb*(j-1)
        else
            thisstart = pos_on_chri[j]
        end
        thisend  = thisstart + window_size_Mb
        push!(window_chr,i)
        push!(window_pos_start,thisstart)
        push!(window_pos_end,thisend)
        snps_window = thisstart .<= pos_on_chri .< thisend
        snps_window_sizej = sum(snps_window)
        push!(window_size_nSNPs,snps_window_sizej)
        if sum(snps_window)!=0
            push!(window_snp_start,pos_on_chri[findfirst(snps_window)])
            push!(window_snp_end,pos_on_chri[findlast(snps_window)])
            push!(window_column_start,index_start)
            push!(window_column_end,index_start+snps_window_sizej-1)
        else #empty windows exist in non-sliding window; no empty window in sliding windows
            push!(window_snp_start,0)
            push!(window_snp_end,0)
        end
        if sliding_window == false
            index_start += snps_window_sizej
        else
            index_start += 1
        end
      end
    end

    out=[]

    if GWAS == true
        if length(threshold) == 1
            println("Compute the posterior probability of association of the genomic window that explains more than ",threshold," of the total genetic variance.")
        elseif threshold == "window specific"
            println("Thresholds for each window are determined based on the window size")
        end
        for i in 1:length(marker_effects_file)
            #using marker effect files
            output            = readdlm(marker_effects_file[i],',',header=true)[1]
            nsamples,nMarkers = size(output)
            nWindows          = length(window_size_nSNPs)
            winVarProps       = zeros(nsamples,nWindows)
            winVar            = zeros(nsamples,nWindows)
            #window_mrk_start ID and window_mrk_end ID are not provided now
            X = (typeof(mme) <: Array ? mme : mme.output_genotypes)
            # calculate window-specific threshold
            if threshold == "window specific"
                threshold0 = window_size_nSNPs ./ nMarkers
                threshold = reshape([threshold0[div(i,nsamples)+1] for i=0:nsamples*nWindows-1], nsamples, nWindows)
            else threshold = threshold
            end

            for i=1:nsamples
                α = output[i,:]
                genVar = var(X*α)
                for winj = 1:length(window_column_start)
                  wStart = window_column_start[winj]
                  wEnd   = window_column_end[winj]
                  BV_winj= X[:,wStart:wEnd]*α[wStart:wEnd]
                  var_winj = var(BV_winj)
                  winVar[i,winj]      = var_winj
                  winVarProps[i,winj] = var_winj/genVar
                end
            end
            winVarProps[isnan.(winVarProps)] .= 0.0 #replace NaN caused by situations no markers are included in the model
            WPPA, prop_genvar = vec(mean(winVarProps .> threshold,dims=1)), vec(mean(winVarProps,dims=1))
            prop_genvar = round.(prop_genvar*100,digits=2)
            winVarmean = vec(mean(winVar,dims=1))
            winVarstd  = vec(std(winVar,dims=1))

            srtIndx = sortperm(WPPA,rev=true)
            outi = DataFrame(trait  = fill(i,length(WPPA))[srtIndx],
                            window = (1:length(WPPA))[srtIndx],
                            chr    = window_chr[srtIndx],
                            wStart = window_pos_start[srtIndx],
                            wEnd   = window_pos_end[srtIndx],
                            start_SNP = window_snp_start[srtIndx],
                            end_SNP   = window_snp_end[srtIndx],
                            numSNP  = window_size_nSNPs[srtIndx],
                            estimateGenVar  = winVarmean[srtIndx],
                            stdGenVar     = winVarstd[srtIndx],
                            prGenVar = prop_genvar[srtIndx],
                            WPPA     = WPPA[srtIndx],
                            PPA_t  = cumsum(WPPA[srtIndx]) ./ (1:length(WPPA)))
             push!(out,outi)
        end
    end
    if genetic_correlation == true && length(marker_effects_file) ==2
            #using marker effect files
            output1           = readdlm(marker_effects_file[1],',',header=true)[1]
            output2           = readdlm(marker_effects_file[2],',',header=true)[1]
            nsamples,nMarkers = size(output1)
            nWindows          = length(window_size_nSNPs)
            gcov              = zeros(nsamples,nWindows)
            gcor              = zeros(nsamples,nWindows)
            #window_mrk_start ID and window_mrk_end ID are not provided now
            X = (typeof(mme) <: Array ? mme : mme.output_genotypes)
            for i=1:nsamples
                α1 = output1[i,:]
                α2 = output2[i,:]
                for winj = 1:length(window_column_start)
                  wStart = window_column_start[winj]
                  wEnd   = window_column_end[winj]
                  Xwinj  = X[:,wStart:wEnd]
                  BV_winj1= Xwinj*α1[wStart:wEnd]
                  BV_winj2= Xwinj*α2[wStart:wEnd]
                  gcov[i,winj] = cov(BV_winj1,BV_winj2)
                  gcor[i,winj] = cor(BV_winj1,BV_winj2)
                end
            end
            gcov[isnan.(gcov)] .= 0.0
            gcor[isnan.(gcor)] .= 0.0
            gcovmean,gcovstd = vec(mean(gcov,dims=1)),vec(std(gcov,dims=1))
            gcormean,gcorstd = vec(mean(gcor,dims=1)),vec(std(gcor,dims=1))
            outi = DataFrame(trait  = fill("cor(t1,t2)",length(gcormean)),
                            window = (1:length(gcormean)),
                            chr    = window_chr,
                            wStart = window_pos_start,
                            wEnd   = window_pos_end,
                            start_SNP = window_snp_start,
                            end_SNP   = window_snp_end,
                            numSNP   = window_size_nSNPs,
                            estimate_cov = gcovmean,
                            std_cov      = gcovstd,
                            estimate_cor = gcormean,
                            std_cor      = gcorstd)
             push!(out,outi)
    end
    return output_winVarProps ? (Tuple(out),winVarProps) : Tuple(out)
end

GWAS2 = GWAS_mod(model,work_path*"controlled_map.csv",result_path*"MCMC_samples_marker_effects_y.txt",
    header=true,window_size="1 Mb",GWAS=true,threshold="window specific")
CSV.write("window_1MB_mod.csv", GWAS2[1])

CSV.write("Pheno_data.csv",Pheno_data)
CSV.write("QTL_data.csv",QTL_data)

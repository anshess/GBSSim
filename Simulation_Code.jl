## This code simulates a population based on provided population parameters and an optional pedigree file.
## Written by Melanie and Andrew Hess

using Distributions
# This file includes all the functions needed for this simulation
include("180309_SimFunc.jl")

# These are the parameters for the population:
# numFounders is the initial number of individuals in the population.
# You can specify the number of chromosomes (numChrs) but at this stage they must all be the same length (chrLength_cM).
numFounders = 500
numChrs = 1
chrLength_cM = 100

# Parameters for the number of SNPs and QTL to simulate. Some will be fixed during the simulation due to sampling/drift,
# therefore we need to simulate more SNPs(numSNP) than the number we want to keep (numKeepSNP) - same for QTL. SNP
# parameters are for each chromosome, while QTL parameters are across all chromosomes.
numSNP = 6000
numKeepSNP = 4000
numQTL = 160
numKeepQTL = 100

# This is the code to simulate the population structure. In our example we simulate 100 founders (numFounders). We then
# call the changingPopSize function whose input parameters are the set of founder animals, the number of animals you want
# to end up with (1000), the number of generations to get to this number (100) and the chromosome length and number of
# chromosomes. Following this we keep population size constant for 900 generations through the constantPopSize function.
founders = generateFounders(numFounders,numChrs);
indCount = [numFounders+1];
inc = changingPopSize(founders,1000,100,chrLength_cM,numChrs);
off = constantPopSize(inc,900,chrLength_cM,numChrs);
inc1 = inc[1:500];
off1 = constantPopSize(inc1,2,chrLength_cM,numChrs);
incX = inc[501:1000];
offX = constantPopSize(incX,5,chrLength_cM,numChrs);
inc2 = offX[1:250];
off2 = constantPopSize(inc2,2,chrLength_cM,numChrs);
incY = offX[251:500];
offY = constantPopSize(incY,5,chrLength_cM,numChrs);
inc3 = offY[1:125];
off3 = constantPopSize(inc3,2,chrLength_cM,numChrs);
inc4 = offY[126:250];
off4 = constantPopSize(inc4,2,chrLength_cM,numChrs);

# This is to simulate from a given pedigree - optional - if you don't use this then run "animals = off" so the code below
#  runs properly
ped = readdlm("smallPed.ped",',',Int64,header=true)[1];
#animals = samplePedigree(ped,off,chrLength_cM,numChrs);
animals = [off1;off2;off3;off4];
#animals = off

# Code to get the SNP and QTL allele frequencies
snpAF = [getAlleleFreq(numSNP) for c in 1:numChrs];
qtlChr = sample([1:numChrs...],numQTL);
numQTLPerChr = [sum(qtlChr .== c) for c in 1:numChrs];
qtlAF = [getAlleleFreq(numQTLPerChr[c]) for c in 1:numChrs];
makeFounderSNPs(founders,snpAF);
makeFounderQTL(founders,qtlAF);

# Code to get the SNP and QTL positions
snpPos = [sampleSNPPosition(size(snpAF[c],1),chrLength_cM,[0.25,0.6,0.15],[0.1,0.8,0.1]) for c in 1:numChrs];
qtlPos = [sample([1:(chrLength_cM*1000000)...],size(qtlAF[c],1)) for c in 1:numChrs];

# This code gets the genotypes for the animals you are interested in - it can take a long time to run depending on your
# population size and number of SNPs/chromosomes
@time fillHaplotypes(animals,founders,numChrs,snpPos,qtlPos)

# This gets the genotypes for your animals of interest for the number of SNPs you specified and writes them out to the
# file - this file name can be changed to whatever you want
snp = getSNPGenotypes(animals);
af = [mean(snp[:,i])/2 for i in 2:size(snp,2)];
goodSNP = find((af .> 0.000) .& (af .<1.000));
keepSNPs = []
for c in 1:numChrs
    startPos = (c-1)*numSNP+1
    endPos = c*numSNP
    goodSNPChr = goodSNP[(goodSNP .>= startPos) .& (goodSNP .<= endPos)]
    keepSNPChr = sort(goodSNPChr[randperm(length(goodSNPChr))[1:numKeepSNP]])
    keepSNPs = [keepSNPs;keepSNPChr]
end
writedlm("dataset1_trueSNPGenos.txt",snp[:,([1;(keepSNPs+1)])])

# This writes out position information for your SNPs
snpPos2 = [0];
snpChr = [0];
for i in 1:numChrs
    snpPos2 = [snpPos2;snpPos[i]]
    snpChr = [snpChr;Array{Int64}(ones(length(snpPos[i]))*i)]
end;
snpPos2 = snpPos2[2:end];
snpChr = snpChr[2:end];
snpDat = [["SNP"*string(i) for i in 1:(numKeepSNP*numChrs)] snpChr[keepSNPs] snpPos2[keepSNPs]];
writedlm("dataset1_snpInfo.txt",snpDat)
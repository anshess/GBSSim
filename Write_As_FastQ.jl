## Code that takes the true SNP genotypes file from the Simulation and outputs a fastq file of "sequenced" reads.
## Written by Melanie Hess.
using Distributions

# Files to read in:
# genos is the file with the true SNP genotypes from the simulation
# seqs is the file of sequences obtained from the getSequences code
# idSeqs is the file with barcodes for each individual - the first individual in the genos file will be allocated to
#    the first barcode in the file etc.
# restrictionEnzyme is the restriction enzyme to use. For now I have just used a single digest of ApeKI
genos = readdlm("dataset1_trueSNPGenos.txt",'\t',Int64,header=false)
seqs = readdlm("sequences.txt",'\t',header=false)
idSeqs = readdlm("IlluminaSeqs.txt",header=false)
restrictionEnzyme = ["CAGC";"G"];

# Mean depth to sequence at - this can be changed
meandepth=6

# Sample the number of reads for each locus for each animal
numAns = size(genos,1)
numSNPs = size(genos,2)-1
snps=rand(Gamma(meandepth*0.5,2),numSNPs)
anis=rand(Gamma(meandepth*3,1/3),numAns)
anis = reshape(anis,numAns,1)
x = (anis * snps')./meandepth
reads = [rand(NegativeBinomial(x[i,j],0.5),1)[1] for i in 1:numAns, j in 1:numSNPs]

# Sample the number of reads for each allele at each locus
nA = [0 for i in 1:numSNPs, j in 1:numAns];
for j in 1:numAns
   genotype = genos[j,2:end]
   numReads = reads[j,:]
   success = Array{Int64}(round.(100*meandepth.*(genotype/2)))
   failure = Array{Int64}(round.(100*meandepth.*(abs.(2-genotype)/2)))
   dist = [Hypergeometric(success[i],failure[i],numReads[i]) for i in 1:size(genotype,1)]
   nA[:,j] = [rand(dist[i],1)[1] for i in 1:size(genotype,1)]
   println("Done animal $j of $numAns")
end
nA = nA'
nB = reads-nA;

# Translate the allele counts into reads that will be sequenced i.e. barcode, restriction site, sequence, restriction site
# NB1 each sequence has restriction cut site at each end - this can be changed if desired
# NB2 the locus is always in the center of the 65bp read - this can also be changed easily
nucl = ["A";"C";"G";"T"]
io = open("tmpOut.txt", "a")
for i in 1:numSNPs
    a1 = rand(nucl)
    a2 = rand(nucl[(nucl .!= a1)])
    s = seqs[i]
    allele1 = restrictionEnzyme[1]*s[1:32]*a1*s[34:end]*restrictionEnzyme[2]
    allele2 = restrictionEnzyme[1]*s[1:32]*a2*s[34:end]*restrictionEnzyme[2]
    a = nA[:,i]
    b = nB[:,i]
    ida = idSeqs.*allele1
    idb = idSeqs.*allele2
    writedlm(io,collect(Iterators.flatten([[repmat([ida[j]],a[j]);repmat([idb[j]],b[j])] for j in 1:numAns])))
end
close(io)

# Reads are shuffled as would be seen in practice and written out in fastq format
# NB the quality score is high for all positions at the moment - can be altered but need to consider the most appropriate
#    way to do this i.e. it may not be desirable to just change the scores without changing any sequence
totalReads = readdlm("tmpOut.txt",header=false)
totalReads = shuffle(totalReads[:,1])
io = open("dataset1_simulatedReads.fastq","w")
for i in 1:size(totalReads,1)
    writedlm(io,["@ Read $i\n"*totalReads[i]*"\n+\n"*repeat("F",length(totalReads[i]))],quotes=false)
end
close(io)
rm("tmpOut.txt")
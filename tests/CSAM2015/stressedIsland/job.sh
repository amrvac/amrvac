#$ -N ff-isl-k4000
#$ -pe ib 64
#$ -l h_rt=48:00:00
#$ -l h_vmem=2000M
#$ -l placement=good
#$ -l h=!h1s1b15n1
#$ -M o.porth@leeds.ac.uk
#$ -cwd -V
#$ -q mhd1.q 
mpirun ./amrvac >>output/out

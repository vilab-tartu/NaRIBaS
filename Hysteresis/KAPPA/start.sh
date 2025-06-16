export OMP_NUM_THREADS=4

conda init
conda activate mds

timestamp=`date +%Y-%m-%d-%H-%M`

nohup ./naribas Kappa_input_STEP1-preparation.sh > ${timestamp}.out
#nohup ./naribas Kappa_input_STEP2-run.sh > ${timestamp}.out

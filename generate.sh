#!/bin/bash

for i in $(seq 0 100); do
cat > "job${i}.sh" <<EOF
#!/bin/bash
#SBATCH -A ccsd
#SBATCH -p burst
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -c 1
#SBATCH -J dqmc
#SBATCH --mem=100G
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cosdis123@gmail.com

cd /home/9pm/ED_qbasis/qbasis/bin/1d/16orb
export OMP_NUM_THREADS=32

time ./chain_Hubbard_tbc_G_calculation.x 8 8 ${i} 0
EOF
chmod +x "job${i}.sh"
done

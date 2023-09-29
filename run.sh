Rscript export-data.R
Rscript aiptasia.R

gcc -o3 optim-csv.c -lm -o optim-csv.ce

./optim-csv.ce cycle-AL.txt 0.1 0 > tmp1 &
./optim-csv.ce cycle-RES.txt 0.1 0 > tmp2 &
./optim-csv.ce data-10d.txt 0.1 0 > tmp3 &
./optim-csv.ce data-140d.txt 0.1 0 > tmp4 &
./optim-csv.ce data-170d.txt 0.1 0 > tmp5 &
./optim-csv.ce data-200d.txt 0.1 0 > tmp6 &
./optim-csv.ce data-21d.txt 0.1 0 > tmp7 &

./optim-csv.ce data-apo.txt 0.1 0 > tmp8 &
./optim-csv.ce data-sym.txt 0.1 0 > tmp9 &

./optim-csv.ce data-170d-170.txt 0.1 0 > tmp10 &
./optim-csv.ce data-170d-200.txt 0.1 0 > tmp11 &

./optim-csv.ce data-21d-sep23.txt 0.1 0 > tmp12 &

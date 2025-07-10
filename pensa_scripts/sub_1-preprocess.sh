NAME=$1

rm preproc_output.txt
rm preproc_error.txt

sbatch -p rondror -t 02:00:00 --mem=80G --wrap='python scripts/1-preprocess.py' -o preproc_output.txt -e preproc_error.txt

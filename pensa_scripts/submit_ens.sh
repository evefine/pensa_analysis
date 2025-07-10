rm job_output.txt
rm job_error.txt

sbatch -p rondror -t 02:00:00 --mem=80G --wrap='python ensemble_comp.py' -o job_output.txt -e job_error.txt


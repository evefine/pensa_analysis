
sbatch -p rondror -t 02:00:00 --mem=80G --wrap="bash scripts/2-compare.sh" -o compare_output.txt -e compare_error.txt

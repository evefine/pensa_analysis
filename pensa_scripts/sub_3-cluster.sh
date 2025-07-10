rm KOR_arr_Gi/cluster_output.txt
rm KOR_arr_Gi/cluster_error.txt

NAME=$1

sbatch -p rondror -t 02:00:00 --mem=80G --wrap="bash 3-cluster.sh $NAME" -o KOR_arr_Gi/cluster_output.txt -e KOR_arr_Gi/cluster_error.txt
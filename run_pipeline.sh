cd /public/home/liunangroup/liangyan/project/bcl11a/20250721.znf_RNA-seq

# TEcount pipeline
thread="1"

for i in $(cat zzsample.txt);do
if [[ -e TEcount/$i.cntTable ]];then
    echo "skipping $i"
else
    sbatch -J job -p cu -c $thread -o %j.out -e %j.err --wrap="
        sh RNA-seq_TE.sh $i $thread
    "
fi
done

# DESeq2 pipeline
thread="16"

sbatch -J job -p cu -c $thread -o %j.out -e %j.err --wrap="
    python DiffExp.py
"
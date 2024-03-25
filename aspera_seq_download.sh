# sudo apt install xlsx2csv
# 下载aspera
wget https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/092u0/0/ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz
# 解压并安装aspera
tar xvf ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.10.0.180973-linux-g2.12-64.sh
# 添加环境变量
echo 'export PATH=~/.aspera/connect/bin/:$PATH' >> ~/.bashrc
# 将metadata_all.xlsx转换为csv格式
xlsx2csv metadata_all.xlsx > metadata_all.csv
# 提取metadata_all.csv的最后一列到aspera.txt
awk -F',' '{print $NF}' metadata_all.csv > aspera.csv
# 删除aspera.csv的第一行
sed -i '1d' aspera.csv
# 将aspera.csv中的;替换为换行符
sed -i 's/;/\n/g' aspera.csv
# 删除aspera.csv中每行的ftp.sra.ebi.ac.uk
sed -i 's/ftp.sra.ebi.ac.uk//g' aspera.csv
# 序列下载
ascp -v -QT -l 400m -P33001 -k1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp --file-list aspera.csv ./seq #用aspera对序列文件进行批量下载
# 对metadata_all.csv只保留第一列和最后一列
awk -F',' '{print $1","$NF}' metadata_all.csv > rename.csv
# 删除rename.csv的第一行
sed -i '1d' rename.csv
# 将rename.csv中的;替换为换行符
sed -i 's/;/\n,/g' rename.csv
# 将rename.csv中第一列的空值替换为前一行的值
awk -F',' '{if($1=="") print a","$2;else print $1","$2; a=$1}' rename.csv > rename2.csv
# 若rename2.csv中第二列以_1.fastq.gz结尾，则第一列末尾加上_1，若以_2.fastq.gz结尾，则第一列末尾加上_2
awk -F',' '{if($2~/_1.fastq.gz$/) print $1"_1,"$2;else if($2~/_2.fastq.gz$/) print $1"_2,"$2;else print $1","$2}' rename2.csv > rename3.csv
# 保留rename3.csv第二列每行最后一个/与.fastq.gz之间的字符
awk -F',' '{split($2, a, "/"); split(a[length(a)], b, ".fastq.gz"); print $1","b[1]".fastq.gz"}' rename3.csv > rename4.csv
# 删除rename4.csv中的.fastq.gz
sed -i 's/.fastq.gz//g' rename4.csv
# 将rename4.csv中的第一列与第二列互换位置
# Swap the first and second columns in rename4.csv
awk -F',' '{print $2","$1}' rename4.csv > rename5.csv
cd seq
# 将rename5.csv换为以制表符分隔的rename5.txt
awk -F',' '{print $1"\t"$2}' ../rename5.csv > ../rename5.txt
# 序列重命名
awk '{system("mv "$1".fastq.gz "$2".fastq.gz")}' ../rename5.txt #序列文件重命名

#!/bin/bash -el
#PBS -l nodes=1:ppn=1,mem=24gb,walltime=72:00:00
#PBS -m abe

cd /home/support/jballer/Seelig

module load python3

python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy2-[20170419-TRIM_1].fastqsanger" -o Library.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy3-[20170419-TRIM_2].fastqsanger" -o ResinLibrary.pep -d -m 8 
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy4-[20170419-TRIM_3].fastqsanger" -o SpeB400R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy5-[20170419-TRIM_4].fastqsanger" -o SpeB400R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy6-[20170419-TRIM_5].fastqsanger" -o SpeB20R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy7-[20170419-TRIM_6].fastqsanger" -o SpeB20R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy8-[20170419-TRIM_7].fastqsanger" -o SpeB1R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy9-[20170419-TRIM_8].fastqsanger" -o SpeB1R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy10-[20170419-TRIM_9].fastqsanger" -o Fxa400R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy11-[20170419-TRIM_10].fastqsanger" -o Fxa400R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy12-[20170419-TRIM_11].fastqsanger" -o Fxa20R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy13-[20170419-TRIM_12].fastqsanger" -o Fxa20R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy14-[20170419-TRIM_13].fastqsanger" -o Fxa1R1.pep -d  -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy15-[20170419-TRIM_14].fastqsanger" -o Fxa1R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy16-[20170419-TRIM_15].fastqsanger" -o ADAM1750R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy17-[20170419-TRIM_16].fastqsanger" -o ADAM1750R2.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy18-[20170419-TRIM_17].fastqsanger" -o ADAM171R1.pep -d -m 8
python3 ProcessingScript.py "/home/support/jballer/Seelig/BoData/Galaxy19-[20170419-TRIM_18].fastqsanger" -o ADAM171R2.pep -d -m 8

#A171R1.fastqsanger   FXa1R1.fastqsanger   FXa400R1.fastqsanger  ResinLibrary.fastqsanger  SpeB1R2.fastqsanger   SpeB20R2.fastqsanger   SpeB400R2.fastqsanger
#A1750R1.fastqsanger  FXa20R1.fastqsanger  Library.fastqsanger   SpeB1R1.fastqsanger       SpeB20R1.fastqsanger  SpeB400R1.fastqsanger


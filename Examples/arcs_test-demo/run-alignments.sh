#!/bin/bash
PATH=/home/rwarren/bin/pigz-2.3.3:$PATH
PATH=/gsc/btl/linuxbrew/bin:$PATH
PATH=/home/lcoombe/bin/samtools-0.1.16:$PATH

pigz --decompress --stdout SOMEBAMFILE | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S+)-/){$flag=1;print "$1_$2\n";}else{$flag=0;}}else{print "$_\n" if($flag);}' > CHROMIUM_interleaved.fastq
bwa mem -t32 YOURDRAFT-renamed.fa -p CHROMIUM_interleaved.fastq | samtools view -Sb - | samtools sort -n - ./CHROMIUM-sorted

#!/bin/bash
/home/rwarren/bin/pigz-2.3.3/unpigz -c SOMEBAMFILE | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print "$_\n" if($flag);}' > CHROMIUM_interleaved.fastq
/gsc/btl/linuxbrew/bin/bwa mem -t32 YOURDRAFT-renamed.fa -p CHROMIUM_interleaved.fastq | /gsc/btl/linuxbrew/bin/samtools view -Sb - | /home/lcoombe/bin/samtools-0.1.16/samtools sort -n - ./CHROMIUM-sorted

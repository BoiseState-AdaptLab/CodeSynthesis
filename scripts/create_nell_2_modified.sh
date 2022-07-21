#!/bin/bash

wget https://s3.us-east-2.amazonaws.com/frostt/frostt_data/nell/nell-2.tns.gz
gzip -d nell-2.tns.gz
echo "# This is a rank 3 test sparse tensor in FROSTT file format," > nell-2-modified.tns
echo "# extended with two meta data lines:" >> nell-2-modified.tns
echo "#   rank nnz" >> nell-2-modified.tns
echo "#   dims (one per rank)" >> nell-2-modified.tns
echo "#" >> nell-2-modified.tns
echo "# see http://frostt.io/tensors/file-formats.html" >> nell-2-modified.tns
echo "#" >> nell-2-modified.tns
echo "# This tensor represents the "B" input to the MTTKRP kernel:" >> nell-2-modified.tns
echo "# http://tensor-compiler.org/docs/data_analytics/index.html" >> nell-2-modified.tns
echo "3 76879419" >> nell-2-modified.tns
echo "12092 9184 28818" >> nell-2-modified.tns
cat nell-2.tns >> nell-2-modified.tns

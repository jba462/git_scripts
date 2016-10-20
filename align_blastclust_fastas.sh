#!/bin/bash

for i in .

    do

        cd ${i}

        for j in *_for_mb

            do

                cd ${j}

                for k in *.fasta

                    do

                        mafft --localpair --anysymbol --maxiterate 1000 --thread 4 ${k} > ${k/.fasta/}.mafft.fas

                    done

                cd ..

            done

        #cd ..

    done

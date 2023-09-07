#!/bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=salmon_merge
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --mail-user=mmarsbrisbin@whoi.edu
#SBATCH --mem=124gb
#SBATCH --time=24:00:00
#SBATCH --output=salmon_merge_%j.log
#SBATCH --error=salmon_merge_%j.err


salmon quantmerge --column NumReads --quants {HM102_quant,HM107_quant,HM115_quant,HM119_quant,HM125_quant,HM86_quant,HM89_quant,HM95_quant,HM103_quant,HM111_quant,HM116_quant,HM120_quant,HM126_quant,HM87_quant,HM92_quant,HM96_quant,HM104_quant,HM112_quant,HM117_quant,HM121_quant,HM127_quant,HM88_quant,HM94_quant,HM97_quant} \
-o salmon.merged.numreads

salmon quantmerge --column TPM --quants {HM102_quant,HM107_quant,HM115_quant,HM119_quant,HM125_quant,HM86_quant,HM89_quant,HM95_quant,HM103_quant,HM111_quant,HM116_quant,HM120_quant,HM126_quant,HM87_quant,HM92_quant,HM96_quant,HM104_quant,HM112_quant,HM117_quant,HM121_quant,HM127_quant,HM88_quant,HM94_quant,HM97_quant} \
-o salmon.merged.TPM

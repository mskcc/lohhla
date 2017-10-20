#!/bin/bash

module load BEDTools/2.26.0-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3
module load novoalign/3.07.00
module load TracerX-Picard-GATK/0.1-Java-1.7.0_80
module load Jellyfish/2.2.6-foss-2016b


Rscript /camp/lab/swantonc/working/rosentr/scripts/lohhla/HLAlossFullScript_updated__user_opt_tmp_RR.R --patientId example --outputDir /camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/test/test-example-file/example-out/ --normalBAMfile /camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/test/test-example-file/bam/example_BS_GL_sorted.bam --BAMDir /camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/test/test-example-file/bam/  --hlaPath /camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/test/test-example-file/hla-types/hlas --HLAfastaLoc /farm/home/lr-tct-lif/wilson52/installs/polysolver/data/abc_complete.fasta --CopyNumLoc /camp/lab/swantonc/working/rosentr/projects/PolySolverLOH/test/test-example-file/ascat/solutions.txt --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir /camp/apps/eb/software/TracerX-Picard-GATK/0.1-Java-1.7.0_80/bin/ --novoDir /camp/apps/eb/software/novoalign/3.07.00/bin/





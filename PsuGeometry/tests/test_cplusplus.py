

import subprocess as sub

dataT = 'x'
pig = sub.Popen(["/home/rachel/Documents/Bioinformatics/BbkProject/Project2/CPlusPlusCode/CMake/CMakeElectronDensity/CMakeElectronDensity",dataT],stdout=sub.PIPE)
result = pig.communicate(input=b"this is sample text.\n")
print(result[0])
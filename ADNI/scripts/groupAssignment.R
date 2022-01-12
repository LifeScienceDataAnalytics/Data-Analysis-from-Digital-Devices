boot_strenwithThreshold <- read_delim("boot_strenwithThreshold.csv", ";", escape_double = FALSE, trim_ws = TRUE)

MMSE = unique(grep("MMSE", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
FAQ = unique(grep("FAQ", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
Imaging =  unique(grep("volume|imaging", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
CSF = unique(grep("csf", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
Amyloid = unique(grep("Amyloid", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
Genetic = unique(grep("APOE4|PHS", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
Diagnosis = unique(grep("DX", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
AR = unique(grep("^AR", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
Motor = unique(grep("^Motor", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to), value = TRUE))
BITDOT = unique(grep("BIT", c(boot_strenwithThreshold$from,boot_strenwithThreshold$to) , value = TRUE))
Demog = unique(grep("AGE|PTMARRY|PTGENDER", c(boot_strenwithThreshold$from, boot_strenwithThreshold$to), value = TRUE))
CognitiveDomains = unique(grep("PerceptualMotorCoordination|ComplexAttention|CognitiveProcessingSpeed|Inhibition|Flexibility|VisualPerception|Planning|ProspectiveMemory|SpatialMemory",c(boot_strenwithThreshold$from, boot_strenwithThreshold$to), value = TRUE))

groupAssignment = cbind.data.frame(c("MMSE", "FAQ", "Imaging", "CSF", "Amyloid", "Genetic",
                                     "Diagnosis", "AR", "Motor", "BITDOT", "Demog", "CognitiveDomains"),
                                   c(paste(MMSE, collapse=", "),paste(FAQ, collapse=", "), paste(Imaging, collapse=", "), paste(CSF, collapse=", "), paste(Amyloid, collapse=", "), 
                                     paste(Genetic, collapse=", "), paste(Diagnosis, collapse=", "), paste(AR, collapse=", "),paste(Motor, collapse=", "), paste(BITDOT, collapse=", "),
                                     paste(Demog, collapse=", "),paste(CognitiveDomains, collapse=", ")))

colnames(groupAssignment) = c("groupNames", "features")

write.csv(groupAssignment, "groupAssignment_ADNI.csv")

# DVSensor – a tool for creating DART VADAR Sensors

DVSensor is a software tool which allows you to create DART VADAR sensors for any mRNA targets. DART VADAR sensors, 
first described by Gayet et al. (Gayet, R.V., Ilia, K., Razavi, S. et al. Autocatalytic base editing for RNA-responsive 
translational control. Nat Commun 14, 1339 (2023). https://doi.org/10.1038/s41467-023-36851-z), are mRNA constructs 
which facilitate conditional mRNA translation based on the 
detection of target mRNA molecules. The DART VADAR mRNA can only be translated after it hybridizes with its target 
mRNA, which activates the sensor and allows the translation of an encoded payload gene. DVSensor takes a target mRNA 
sequence as input, either in FASTA or GenBank format, and generates possible DART VADAR sensor sequences as output. 
The sensor sequences can then be cloned into a suitable DART VADAR expression vector, such as our pASTERISK. The 
software offers a variety of settings for controlling how the sensors are generated. In addition, it comes with a 
built-in feature for evaluating the specificity of the generated sensor sequences. By querying mRNA databases using 
BLAST, it can identify potential off-target mRNA transcripts other than the intended target that may also be able to 
activate the sensor. A menu for documentation and help is available in the software as well. DVSensor runs locally as 
a web app inside the browser and is available for both Linux and Windows. Alternatively, it could also be hosted on 
a web server.  

This README only gives a brief overview over the installation and usage. For a thorough documentation, visit the
[Software](https://2023.igem.wiki/bielefeld-cebitec/software) page of the iGEM Bielefeld-CeBiTec 2023 Team.



## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might
be unfamiliar with (for example your team wiki). A list of Features or a Background subsection can also be added here.
If there are alternatives to your project, this is a good place to list differentiating factors.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew.
However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing
specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a
specific context like a particular programming language version or operating system or has dependencies that have to be
installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of
usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably
include in the README.

## Contributing
This software will not be maintained on this GitLab repository beyond the duration of the iGEM 2023 competition.  
A clone of this repository can be found on github (...), which may or may not be maintained in the future.  
Since this software is licensed under the MIT license, you are free to copy, modify, and distribute
it without restriction. 

## Authors and acknowledgment
The software was developed by Daniel Prib (contact: dpx0@mailbox.org, github: https://github.com/dpx0).  
The following image assets used in the software were obtained from [flaticon.com](https://www.flaticon.com/)
under a free-to-use license:  
[Arrow icon created by Kirill Kazachek - Flaticon](https://www.flaticon.com/free-icons/arrow)  
[Book icon created by Good Ware – Flaticon](https://www.flaticon.com/free-icons/book)  
[Duplicate icon created by Phoenix Group – Flaticon](https://www.flaticon.com/free-icons/duplicate)  
[Document icon created by Freepik – Flaticon](https://www.flaticon.com/free-icons/document)  
[Home button icon created by Freepik – Flaticon](https://www.flaticon.com/free-icons/home-button)  
[Rna icon created by Freepik - Flaticon](https://www.flaticon.com/free-icons/rna)  
[Rna icon created by Rakib Hassan Rahim - Flaticon](https://www.flaticon.com/free-icons/rna)  


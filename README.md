# DVSensor – a tool for creating DART VADAR Sensors

This software tool was created by the Bielefeld-CeBiTec team for the [2023 iGEM competition](https://igem.org/).

## Description
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

This README only gives a brief summary over the installation and usage. For a thorough documentation including
a troubleshooting guideline, visit the
[Software](https://2023.igem.wiki/bielefeld-cebitec/software) page on the iGEM Bielefeld-CeBiTec 2023 team wiki.

## Installation
To use this application, the following software is required:
* Python (version 3.11 or higher)
* Python package: nicegui (version 1.3.9 or higher)
* Python package: Biopython (version 1.81 or higher)  

To install python for your operating system, follow the guides on the official 
[Python website](https://www.python.org/downloads/). Python comes with a 
package manager called pip, which you can use to install python packages. Run the following commands in a 
terminal / command line (cmd.exe on Windows) to install the required packages:  
	* pip install nicegui
	* pip install biopython  

On Linux systems, installing python packages using pip may not work, and requires either an installation with the 
systems package manager, or the use of virtual environments. The official
[Python website](https://packaging.python.org/en/latest/tutorials/installing-packages/) provides 
more information. Also, 
[this thread](https://askubuntu.com/questions/1465218/pip-error-on-ubuntu-externally-managed-environment-×-this-environment-is-extern) 
on the askubuntu forum may be helpful.  

The source code of the application can be found inside the dvsensor/src directory. Download the src directory and 
launch the application by running the main.py file inside a terminal / command line: "python main.py"  

Additionally, there is a build.py file which you can use to create bundled executables with PyInstaller. 
These bundled executables can then be run on any machine without the need for installing python or any packages. 
In order to use this file, install the latest version of the PyInstaller package and run the file in the 
terminal / command line. Unfortunately, we could not provide pre-bundled executables in this repository, 
since the project storage space is limited. If you decide to create a bundled version, you can simply start the 
app by launching the dvsensor.exe file on windows or the dvsensor binary on Linux that is created when 
you run the build.py file.  


## Usage
1. Start the application (...)

2. Open a web browser of your choice and enter this url: http://localhost:8080. You should now see the start page of 
the web app (see image below). **Note:** closing or reloading the browser window will shut down the application. This is 
intended behavior, but requires a restart of the application. You should only use the buttons inside the application 
window to navigate between pages, and not the back/forward buttons of your browser.

3. You are presented with two large buttons which allow you to generate new sensors or read the documentation, 
respectively. Click on the right button to open the documentation, or the left button to upload an mRNA target
sequence.  

![Main Menu](images/01.png "Main Menu")  

4. Upload your target mRNA sequence, either as a FASTA or GenBank record. You can upload a file or paste the
sequence record manually. A maximum record size of 500MB can be uploaded.  

![Sequence Upload](images/02.png "Sequence Upload")  

To upload a file, click on the plus button in the right corner. Select a file and click the upload button. 
To manually enter a sequence, paste the sequence into the text field and click "continue".

![File Upload](images/03.png "File Upload")  

5. The sequence name and NCBI accession number are extracted from the uploaded sequence record. Usually, the sequence 
name gets interpreted as the accession number, which may not be desired. Therefore, you can change both the sequence 
name and accession number. Changing the name has no effect on the computation or the results, but the 
accession number has to be correct in order for BLAST queries to work properly. Click "OK" when you are done.  

![Sequence Info](images/04.png "Sequence Info")  

6. The next menu provides different options for controlling how the sensors are generated and for controlling the BLAST 
queries. The different options are explained on the [team wiki](https://2023.igem.wiki/bielefeld-cebitec/software). 
Click "Run Analysis" when the desired options have been set.  

![Options](images/05.png "Options")  

7. The application will locate the selected triplets inside the selected transcript regions, and generate 123 bp 
sensors centered on those triplets. If the "run blast" option was selected, the trigger sequences (the segments of 
the target mRNA that are centered around a given triplet) will be queried against the specified database using BLAST 
to identify potential mRNA transcripts that may produce an off-target sensor activation. The results are reported 
in a table. The content of the output table is explained in more detail on the 
[team wiki](https://2023.igem.wiki/bielefeld-cebitec/software).  

![Results](images/06.png "Results")  

In the upper left corner, the target mRNA is depicted as a simplified ideogram, with the sizes of the different 
regions corresponding to their relative lengths. Selecting an entry in the table displays a red window in the 
ideogram, which indicates the region in the mRNA that is targeted by that particular sensor.  

![mRNA Ideogram](images/07.png "mRNA Ideogram")  

The "cancel" button in the upper right corner allows you to cancel the running analysis. When the analysis is finished 
or if it was cancelled, a new button appears, which allows you to download the output table as a CSV file.  

![Export table](images/08.png "Export table")  

Clicking the "home" button in the lower left corner brings you back to the start page.

## Contributing
Since this software was developed for the 2023 iGEM competition, it will not be maintained on this GitLab repository 
beyond the duration of the competition. A clone of this repository can be found on github (https://github.com/dpx0) 
where it may be maintained in the future. This software is licensed under the 
[MIT license](https://opensource.org/license/MIT/), so you are free to copy, modify, and distribute it 
without restriction. 

## Authors and acknowledgment
This software was developed by Daniel Prib (contact: dpx0@mailbox.org, github: https://github.com/dpx0).  
The following image assets used in the software were obtained from [flaticon.com](https://www.flaticon.com/)
under a free-to-use license:  
[Arrow icon created by Kirill Kazachek - Flaticon](https://www.flaticon.com/free-icons/arrow)  
[Book icon created by Good Ware – Flaticon](https://www.flaticon.com/free-icons/book)  
[Duplicate icon created by Phoenix Group – Flaticon](https://www.flaticon.com/free-icons/duplicate)  
[Document icon created by Freepik – Flaticon](https://www.flaticon.com/free-icons/document)  
[Home button icon created by Freepik – Flaticon](https://www.flaticon.com/free-icons/home-button)  
[Rna icon created by Freepik - Flaticon](https://www.flaticon.com/free-icons/rna)  
[Rna icon created by Rakib Hassan Rahim - Flaticon](https://www.flaticon.com/free-icons/rna)  


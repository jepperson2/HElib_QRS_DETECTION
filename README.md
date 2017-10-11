# Homomorphic QRS Detection **(hqrsd)**

*This project is being developed as part of my Master's thesis at the University of Malaya*

*Although this project is still in development, feel free to contact me at epperson . jessee at gmail . com*


## 1. What is it, in *one* line?
It is an implementation of a QRS detection algorithm using fully homomorphic encryption (HElib + hbc).


## 2. What is in there?
- It is written in C++ and is cross-platform with Vagrant
- The binary circuit API (as developed by qdm12) in *src/he.cpp*
- The code that is run as an **example** is in *src/main.cpp*
- Some other classes are in *src/helper_functions.cpp*
- The QRS detection algorithm is implemented in *src/QRS_DETECTION.cpp*
- All the other src files are **unit tests** and timing tests for the homomorphic binary operations implemented in *src/he.cpp*
- There is a **Vagrantfile** to setup eveything for you, cross-platform
- There is a **makefile** to build *hqrsd* or setup almost everything for you (depending on your OS).
- There is this complete, detailed and updated **README.md** file.


## 3. What does it run ?
- It runs the code in *main.cpp* which executes the QRS detection algorithm on ECG data found in a file specified in *main.cpp*.
- You can change main.cpp with your code to suit your needs/inputs.


## 4. What does it require ?
- Practically:
    - A Linux/Windows/OS X computer
    - At least 3GB of RAM and 2 CPU cores
    - An internet connection
    - CPU with Hardware virtualization tech ideally (you probably have it don't worry)
- In terms of software (although this is automatically installed), here are the dependencies:


   | Program or Library | Requirement 1 | Requirement 2 | Requirement 3 | Requirement 4 | Requirement 5 |
   | ------------------ | ------------- | ------------- | ------------- | ------------- | ------------- |
   | hqrsd		| g++		| make		| HElib		| c++11		| 		|
   | HElib              | g++           | make          | git           | libboost      | NTL 10.5.0    |
   | NTL 10.5.0         | g++           | make          | GMP 6.1.2     |		|		|
   | GMP 6.1.2          | g++           | make          | m4            | perl          |		|


## 5. Documentation
- This readme file
- Comments in the source code, especially in _QRS_DETECTION.cpp_
- A link to my written thesis to come, once I finish writing it :) 


## 6. Abstract (Draft) ##
Cardiac arrhythmias are associated with heart diseases and are often asymptomatic. Correct detection of arrhythmias helps diagnose and treat heart diseases. The most common method of detecting arrhythmias is electrocardiogram (ECG) monitoring. Some arrhythmias, such as atrial fibrillation, are paroxysmal in nature and cannot be easily monitored in a hospital setting. Healthcare Organizations (HCOs) must be careful, however, when collecting data outside their confines, as the privacy of medical information could be compromised during the transfer, storage, or computation of said data. Secure transfer and storage of ambulatory ECG data is trivially accomplished by use of well-established encryption and transfer protocols, but privacy-preserving computation is non-trivial. Existing solutions for analysis of ECG data require either delayed analysis once apparatus is brought back to the HCO or significant infrastructure to transfer and analyze the data in a secure environment. Using public cloud computing would significantly reduce the amount of infrastructure required for an HCO to analyze ECG data, but allowing third parties unencrypted access to medical data goes against HIPAA regulations and violates patient privacy. Private clouds cannot be feasibly instantiated by small HCOs due to the high set-up and maintenance costs, though they would provide the security guarantees required for privacy-preserving computation. One privacy-preserving solution that allows public cloud servers to compute sensitive medical data is to use fully homomorphic encryption (FHE), which enables computations to be performed on encrypted data without decryption. Using FHE would allow HCOs to transfer, store, and process medical information in untrusted, public clouds while staying within the privacy guidelines of HIPAA because the cloud servers would never see the unencrypted data. This work uses HElib and an API called hbc to create an FHE-based implementation of a QRS complex detection algorithm as an example of ECG data computation.

## 7. How do I run it?

### 7.1 Using Vagrant (easiest, compatible with all, most flexible)
1. Install git on your computer
    - `apt-get install -y git` for Linux machines
    - or download it from [git-scm.com/downloads](https://git-scm.com/downloads)
2. On Windows, have an ssh client or add the **ssh.exe** of `C:\Program Files\Git\usr\bin` to your environment path
2. Install Virtual Box from [virtualbox.org/wiki/Downloads](https://www.virtualbox.org/wiki/Downloads)
3. Install Vagrant from [vagrantup.com/downloads.html](https://www.vagrantup.com/downloads.html)
4. Open a terminal and enter `git clone https://github.com/jepperson2/hqrsd.git`
5. Go to the hqrsd directory with `cd hqrsd`
6. Enter `vagrant up` to launch the virtual machine which will setup and build everything for you. 
This takes about 30 minutes the first time, depending on your connection speed and CPU.
This basically launches an Ubuntu-based virtual machine with only what is necessary for this project.
**WARNING:** If you do not have hardware virtualization, you can still run it but you have to change *trusty64*
 to *trusty32* and *vb.cpus = 2* to *vb.cpus = 1*.
7. Once vagrant up has completed, enter `vagrant ssh` to log in the virtual machine.
8. The working directory *hqrsd* on your host machine is shared with the virtual machine at `/vagrant`.
9. In the virtual machine, enter `cd /vagrant`.
10. What's nice then:
    - You can modify the files on your host machine (like Windows etc.)
    - Changes you make are automatically reflected in the Ubuntu-based virtual machine.
    - Compile hqrsd again with `make hqrsd` in the virtual machine.
    - Run hqrsd with ./hqrsd from the virtual machine or your host machine.
    - **Note:** You can use `make hqrsdNrun` to build and automatically run the main.cpp code.
11. When you are done:
    - Enter `exit` in the virtual machine, bringing you back to your host machine.
    - Enter `vagrant halt` to shutdown the machine. Or enter `vagrant destroy` to delete the machine.
12. To log back in, enter `vagrant up` and it should take about 30 seconds ! (except if you destroy the machine)

### 7.2 Using more manual ways, which don't work for all OSes
1. Make sure you have installed **make**
2. Open a terminal as **root** or **administrator** ideally
3. Setup the necessary libraries
    - With the Makefile provided (only works for **Debian** and **Ubuntu**)
        1. Note: *git, g++, m4, perl, gmp and ntl* will be installed automatically*.
        2. Enter `make HElib` in a terminal in the *hqrsd* directory.
    - Manually (if Vagrant and Makefile are not good for you)
        - Mac OSX
            1. Install Xcode manually or with `xcode-select --install`
            2. Install brew with `ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
            3. Install libraries with `brew install wget git g++ m4 perl libboost`
            4. Download GMP with `curl https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2 > gmp-6.1.2.tar.bz2`
            5. Extract it and go to its directory with `tar -xvjf gmp-6.1.2.tar.bz2 && cd gmp-6.1.2`
            6. Configure it with `./configure`
            7. Build it with `make`
            8. Install it with `make install`
            9. *Optionally*, check it with `make check`
            10. Go back and remove used files with `cd .. && rm -fr gmp-6.1.2*`
            10. Download NTL with `curl http://www.shoup.net/ntl/ntl-10.5.0.tar.gz > ntl-10.5.0.tar.gz`
            11. Extract it and go to its directory with `tar -xvzf ntl-10.5.0.tar.gz && cd ntl-10.5.0/src`
            12. Configure it with `./configure NTL_GMP_LIP=on`
            13. Build it with `make`
            14. Install it with `make install`
            15. Go back and remove used files with `cd ../.. && rm -fr ntl-10.5.0*`
            16. Clone HElib with with `git clone https://github.com/shaih/HElib.git`
            17. Go to its src directory `cd HElib/src`
            18. Build it with `make`
            19. *Optionally*, check it with `make check` and test it with `make test`.
            20. Go back to the working directory with `cd ../..`
        - Other Linux OSes
            1. Install the libaries with (add `sudo` maybe) `apt-get install git g++ m4 perl libboost-all-dev`
            2. Download GMP with `wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2`
            3. Extract it and go to its directory with `tar -xvjf gmp-6.1.2.tar.bz2 && cd gmp-6.1.2`
            4. Configure it with `./configure`
            5. Build it with `make`
            6. Install it with `make install`
            7. *Optionally*, check it with `make check`
            8. Go back and remove used files with `cd .. && rm -fr gmp-6.1.2*`
            9. Download NTL with `wget http://www.shoup.net/ntl/ntl-10.5.0.tar.gz`
            10. Extract it and go to its directory with `tar -xvzf ntl-10.5.0.tar.gz && cd ntl-10.5.0/src`
            11. Configure it with `./configure NTL_GMP_LIP=on`
            12. Build it with `make`
            13. Install it with `make install`
            14. Go back and remove used files with `cd ../.. && rm -fr ntl-10.5.0*`
            15. Clone HElib with with `git clone https://github.com/shaih/HElib.git`
            16. Go to its src directory `cd HElib/src`
            17. Build it with `make`
            18. *Optionally*, check it with `make check` and test it with `make test`.
            19. Go back to the working directory with `cd ../..`
        - Cygwin 32 bit and 64 bit
			- It will fail when you try to install NTL with the NTL_CMP_LIP=on because cygwin does not find the -lgmp libary for some reason.
			- So just switch to use Vagrant. I might be missing something but there is no point digging to deep here I believe.
			- Here would be the procedure:
               1. Close any Cygwin processes running
               2. Download the right Cygwin installer:
                  - [Cygwin 32 bit](http://www.cygwin.com/setup-x86.exe)
               	  - [Cygwin 64 bit](http://www.cygwin.com/setup-x86_64.exe)
               3. Run the Cygwin installer previously downloaded
               4. Click *Next >*, *Next >*, *Next >*, *Next >*, *Next >*, *Next >*
               5. With the help of the search bar, select the following packages (only the Devel is necessary):
                  - git: Distributed version control system
                  - gcc-g++: GNU Compiler Collection (C++)
                  - make: the GNU version of the 'make' utility
                  - m4: GNU implementation of the traditional Unix macro processor
                  - perl: Perl Programming language interpreter
                  - libboost-devel: Boost C++ libraries
                  - gmp: Library for arbitrary precision arithmetic
               6. Click on *Next >*, *Next >* and wait for the installation to finish and then close the window.
               7. Lauch a Cygwin terminal
               8. Download NTL with `wget http://www.shoup.net/ntl/ntl-10.5.0.tar.gz`
               9. Extract it and go to its directory with `tar -xvzf ntl-10.5.0.tar.gz && cd ntl-10.5.0/src`
               10. Configure it with `./configure NTL_GMP_LIP=on`			
               11. Build it with `make`
               12. Install it with `make install`
               13. Go back and remove used files with `cd ../.. && rm -fr ntl-10.5.0*`
               14. Clone HElib with with `git clone https://github.com/shaih/HElib.git`
               15. Go to its src directory `cd HElib/src`
               16. Build it with `make`
               17. *Optionally*, check it with `make check` and test it with `make test`.
               18. Go back to the working directory with `cd ../..`
4. Build hqrsd
    - With the Makefile provided (compatible will **all** platforms).
        1. Build it with `make hqrsd`
    - Manually
        1. Create the directory objects `mkdir -p objects`
        2. Compile the API `g++ -c src/he.cpp -I HElib/src -o objects/he.o`
        3. Compile the helper functions `g++ -c src/helper_functions.cpp -o objects/helper_functions.o`
        4. Compile the various tests
            - `g++ -c src/TEST_GATES.cpp -I HElib/src -o objects/test_gates.o`
            - `g++ -c src/TEST_CIRC_COMB.cpp -I HElib/src -o objects/test_circ_comb.o`
            - `g++ -c src/TEST_CIRC_SEQ.cpp -I HElib/src -o objects/test_circ_seq.o`
            - `g++ -c src/TEST_CIRC_ARITHM.cpp -I HElib/src -o objects/test_circ_arithm.o`
	    - `g++ -c src/QRS_DETECTION.cpp -I HElib/src -o objects/QRS_DETECTION.o
        5. Compile the main.cpp file `g++ -c src/main.cpp -I HElib/src -o objects/main.o`
        6. Compile the objects into *hqrsd* `g++ objects/*.o HElib/src/fhe.a -o hqrsd -L/usr/local/lib -lntl -lgmp -lm`
5. Run hqrsd
    - Run it with `./hqrsd` (Careful about having enough **RAM**)
    - You can also build it and run the new build with `make hqrsdNrun`

	
## 8. RAM considerations IMPORTANT
- To run the default hqrsd program, at least 3GB of RAM is recommended.
- Note that applications which use simpler circuits will require less RAM, as simplier circuits like multiplication only require about 0.7 - 1GB of RAM, whereas the averages circuit takes about 3GB
- For **Vagrant**, you can modify the amount of RAM in the **vb.memory** field, 
  which is set to **7000MB** by default. To monitor the RAM usage, open a new 
  host terminal, go to the working directory and use `vagrant ssh -c htop`.

	  
## 9. CPU considerations for Vagrant
- By default, the Vagrant VM uses 2 cores of your CPU (vb.cpus = 2) so that
  you can run hqrsd and also monitor the RAM with another `vagrant ssh`.
- You can also run more instances of hqrsd if you have more than two cores available.
  With Vagrant, just set vb.cpus to 3 for example, log in with `vagrant ssh` on different
  host terminals and run hqrsd (provided you have enough RAM to run both obviously).

	  
## 10. Remove and uninstall ##

### 10.1 With Vagrant
Just enter `vagrant destroy` from your host machine in the working directory.

### 10.2 Otherwise
Use the makefile and run `make deepclean` which uninstalls and deletes:
- hqrsd
- HElib, NTL, GMP
- perl, m4, git, gcc-g++ and libboost-all-dev and purge them.
Only the makefile will remain in the folder.


## 11. Acknowledgements ##
Credits to **Shai Halevi** for developing and maintaining HElib

Credits to **Quentin McGaw** for developing and maintaining hbc, and for answering my questions about using his API

Thanks to **Prof. Miss Laiha** (University Malaya) for her supervision over my research

Thanks to my beloved **Amy** for her continuous love and support throughout my adventures in Southeast Asia

